# Copyright 2024 Department of Computational Biology for Infection Research - Helmholtz Centre for Infection Research
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import pandas as pd
import sys
import logging
import traceback
import os
import mimetypes
import gzip
import io
import tarfile
import zipfile
from multiprocessing.pool import ThreadPool
from collections import defaultdict
from collections import OrderedDict
from cami_amber import binning_classes
from cami_amber.utils import labels as utils_labels

try:
    import load_ncbi_taxinfo
except ImportError:
    sys.path.append(os.path.dirname(__file__))
    try:
        import load_ncbi_taxinfo
    finally:
        sys.path.remove(os.path.dirname(__file__))


def open_generic(file):
    file_type, file_encoding = mimetypes.guess_type(file)

    try:
        if file_encoding == 'gzip':
            if file_type == 'application/x-tar':  # .tar.gz
                tar = tarfile.open(file, 'r:gz')
                f = tar.extractfile(tar.getmembers()[0])
                return io.TextIOWrapper(f)
            else:  # .gz
                return gzip.open(file, 'rt')
        if file_type == 'application/zip':  # .zip
            f = zipfile.ZipFile(file, 'r')
            f = f.open(f.namelist()[0])
            return io.TextIOWrapper(f)
        else:
            return open(file, 'rt')
    except FileNotFoundError as e:
        logging.getLogger('amber').error('File not found: %s' % file)
        raise e


def open_coverages(file_path):
    if not file_path:
        return pd.DataFrame()
    logging.getLogger('amber').info('Loading coverage file')
    coverages_pd = pd.DataFrame()
    try:
        samples_metadata = read_metadata(file_path)
        for metadata in samples_metadata:
            nrows = metadata[1] - metadata[0] + 1
            df = pd.read_csv(file_path, sep='\t', comment='#', skiprows=metadata[0], nrows=nrows, header=None)
            df.rename(columns=metadata[3], inplace=True)
            df = df[['GENOMEID', 'COVERAGE']]
            df['SAMPLEID'] = metadata[2]['SAMPLEID']
            coverages_pd = pd.concat([coverages_pd, df], ignore_index=True)
    except BaseException as e:
        traceback.print_exc()
        logging.getLogger('amber').critical("File {} not found or malformed. {}".format(file_path, e))
        exit(1)
    return coverages_pd


def load_unique_common(unique_common_file_path, args_keyword):
    if not unique_common_file_path:
        return None
    logging.getLogger('amber').info('Loading list of genomes to be removed')
    with open(unique_common_file_path) as f:
        line = f.readline()
        line_split = line.split('\t')
        f.seek(0)
        if len(line_split) == 1:
            return [line.strip('\n') for line in f]
        elif len(line_split) > 1:
            genome_to_unique_common = {}
            for line in f:
                line_split = line.split('\t')
                genome_to_unique_common[line_split[0]] = line_split[1].strip('\n')
            if args_keyword:
                return [genome_id for genome_id, keyword in genome_to_unique_common.items() if keyword == args_keyword]
            else:
                return list(genome_to_unique_common.keys())
    return None


def load_ncbi_info(ncbi_dir):
    if ncbi_dir:
        logging.getLogger('amber').info('Loading NCBI taxonomy from %s' % ncbi_dir)
        try:
            if os.path.isfile(ncbi_dir):
                taxonomy_df = pd.read_feather(ncbi_dir).set_index('TAXID')
            else:
                taxonomy_df = pd.read_feather(os.path.join(ncbi_dir, 'nodes.amber.ft')).set_index('TAXID')
        except BaseException:
            traceback.print_exc()
            logging.getLogger('amber').info('Preprocessed NCBI taxonomy file not found. Creating file {}'.format(os.path.join(ncbi_dir, 'nodes.amber.ft')))
            taxonomy_df = load_ncbi_taxinfo.preprocess_ncbi_tax(ncbi_dir)
        taxonomy_df = taxonomy_df.astype(dtype={rank: pd.UInt32Dtype() for rank in load_ncbi_taxinfo.RANKS})
        return taxonomy_df


def read_metadata(path_label_tuple):
    if type(path_label_tuple) == tuple:
        file_path_query, label = path_label_tuple

        header = {}
        columns_list = []
        got_column_indices = False

        samples_metadata = []

        with open_generic(file_path_query) as f:
            for i, line in enumerate(f):
                line = line.strip()

                if len(line) == 0 or line.startswith("#"):
                    if got_column_indices and not reading_data:
                        data_start = i + 1
                    continue

                # parse header with column indices
                if line.startswith('@@'):
                    for index, column_name in enumerate(line[2:].split('\t')):
                        columns_list.append(column_name)
                    got_column_indices = True
                    reading_data = False
                    data_start = i + 1

                elif line.startswith('@'):
                    logging.getLogger('amber').info('Found {} in {}'.format(line, label))
                    # parse header with metadata
                    if got_column_indices:
                        samples_metadata.append((data_start, data_end, header, columns_list, file_path_query, label))
                        header = {}
                        columns_list = []
                    key, value = line[1:].split(':', 1)
                    header[key.upper()] = value.strip()
                    got_column_indices = False
                    reading_data = False

                else:
                    reading_data = True
                    data_end = i
        try:
            samples_metadata.append((data_start, data_end, header, columns_list, file_path_query, label))
        except UnboundLocalError as e:
            logging.getLogger('amber').critical("File {} is malformed.".format(file_path_query))
            raise e
        return samples_metadata
    else:

        header = {}
        index_to_column_name = {}
        got_column_indices = False

        samples_metadata = []

        with open(path_label_tuple, 'r') as f:
            for i, line in enumerate(f):
                line = line.strip()

                if len(line) == 0 or line.startswith("#"):
                    if got_column_indices and not reading_data:
                        data_start = i + 1
                    continue

                # parse header with column indices
                if line.startswith('@@'):
                    for index, column_name in enumerate(line[2:].split('\t')):
                        index_to_column_name[index] = column_name
                    got_column_indices = True
                    reading_data = False
                    data_start = i + 1

                elif line.startswith('@'):
                    logging.getLogger('amber').info('Found {}'.format(line))
                    # parse header with metadata
                    if got_column_indices:
                        samples_metadata.append((data_start, data_end, header, index_to_column_name))
                        header = {}
                        index_to_column_name = {}
                    key, value = line[1:].split(':', 1)
                    header[key.upper()] = value.strip()
                    got_column_indices = False
                    reading_data = False

                else:
                    reading_data = True
                    data_end = i
        samples_metadata.append((data_start, data_end, header, index_to_column_name))
        return samples_metadata


def load_sample(metadata):
    columns = ['SEQUENCEID', 'BINID', 'TAXID', 'LENGTH', '_LENGTH']
    logging.getLogger('amber').info('Loading %s of %s' % (metadata[2]['SAMPLEID'], metadata[5]))
    nrows = metadata[1] - metadata[0] + 1
    usecols = [v for v in metadata[3] if v in columns]
    df = pd.read_csv(metadata[4], sep='\t', comment='#', skiprows=metadata[0], nrows=nrows, header=None,
                     names=metadata[3],
                     usecols=usecols,
                     dtype={'SEQUENCEID': pd.StringDtype(), 'BINID': pd.StringDtype(), 'TAXID': pd.UInt32Dtype(),
                            'LENGTH': pd.UInt32Dtype(), '_LENGTH': pd.UInt32Dtype()})
    df.rename(columns={'_LENGTH': 'LENGTH'}, inplace=True)
    return df


def load_binnings(samples_metadata, file_path_query):
    columns = ['SEQUENCEID', 'BINID', 'TAXID', 'LENGTH', '_LENGTH']
    sample_id_to_query_df = OrderedDict()
    for metadata in samples_metadata:
        logging.getLogger('amber').info('Loading ' + metadata[2]['SAMPLEID'])
        nrows = metadata[1] - metadata[0] + 1
        col_indices = [k for k, v in metadata[3].items() if v in columns]
        df = pd.read_csv(file_path_query, sep='\t', comment='#', skiprows=metadata[0], nrows=nrows, header=None, usecols=col_indices)
        df.rename(columns=metadata[3], inplace=True)
        df.rename(columns={'_LENGTH': 'LENGTH'}, inplace=True)
        sample_id_to_query_df[metadata[2]['SAMPLEID']] = df
    return sample_id_to_query_df


def open_query(file_path_query):
    try:
        samples_metadata = read_metadata(file_path_query)
        sample_id_to_query_df = load_binnings(samples_metadata, file_path_query)
    except BaseException as e:
        traceback.print_exc()
        logging.getLogger('amber').critical("File {} not found or malformed. {}".format(file_path_query, e))
        exit(1)
    return sample_id_to_query_df


def get_sample_id_to_num_genomes(sample_id_to_queries_list):
    sample_id_to_num_genomes = {}
    for sample_id in sample_id_to_queries_list:
        if isinstance(sample_id_to_queries_list[sample_id][0], binning_classes.GenomeQuery):
            gs_df = sample_id_to_queries_list[sample_id][0].df
            sample_id_to_num_genomes[sample_id] = gs_df['BINID'].nunique()
        else:
            sample_id_to_num_genomes[sample_id] = 0
    return sample_id_to_num_genomes


def get_rank_to_df(query_df, taxonomy_df, label, is_gs=False):
    if is_gs:
        if 'LENGTH' not in query_df.columns:
            logging.getLogger('amber').critical("Sequences length could not be determined. Please add column LENGTH to gold standard.")
            exit(1)
        cols = ['TAXID', 'LENGTH']
    else:
        cols = ['TAXID']

    rank_to_df = dict()
    logging.getLogger('amber').info('Setting taxa ranks')
    query_df = pd.merge(query_df, taxonomy_df.reset_index(), on=['TAXID']).set_index('SEQUENCEID')

    for rank in load_ncbi_taxinfo.RANKS:
        logging.getLogger('amber').info('Deriving bins of {} at rank: {}'.format(label, rank))
        if query_df[rank].isnull().all():
            continue
        rank_to_df[rank] = query_df[query_df[rank].notnull()][cols + [rank]]
        rank_to_df[rank]['TAXID'] = rank_to_df[rank][rank]
        rank_to_df[rank] = rank_to_df[rank].drop(columns=rank).reset_index()

    return rank_to_df


def load_queries(gold_standard_file, bin_files, labels, options=None, options_gs=None):
    max_workers = min(len(labels) + 1, os.cpu_count() or 1)
    pool = ThreadPool(max_workers)

    try:
        metadata_all = pool.map(read_metadata, zip([gold_standard_file] + bin_files, [utils_labels.GS] + labels))
        pool.close()
    except BaseException:
        logging.getLogger('amber').error('An error occurred. Exiting.')
        exit(1)
    samples_metadata_gs = metadata_all[0]
    samples_metadata_queries = metadata_all[1:]

    if not options:
        options = binning_classes.Options()
    if not options_gs:
        options_gs = binning_classes.Options()
    taxonomy_df = load_ncbi_info(options.ncbi_dir) if options.ncbi_dir else pd.DataFrame()

    sample_id_to_g_queries_list = defaultdict(list)
    sample_id_to_t_queries_list = defaultdict(list)
    sample_id_to_g_gs = defaultdict(list)
    sample_id_to_t_gs = defaultdict(list)

    for metadata in samples_metadata_gs:
        sample_id = metadata[2]['SAMPLEID']
        columns = metadata[3]
        if 'BINID' in columns:
            g_query_gs = binning_classes.GenomeQuery(utils_labels.GS, sample_id, options_gs, metadata, True)
            g_query_gs.gold_standard = g_query_gs
            sample_id_to_g_gs[sample_id] = g_query_gs
            if not options.skip_gs:
                sample_id_to_g_queries_list[sample_id].append(g_query_gs)
        if 'TAXID' in columns and not taxonomy_df.empty:
            t_query_gs = binning_classes.TaxonomicQuery(utils_labels.GS, sample_id, options_gs, metadata, taxonomy_df, True)
            t_query_gs.gold_standard = t_query_gs
            sample_id_to_t_gs[sample_id] = t_query_gs
            if not options.skip_gs:
                sample_id_to_t_queries_list[sample_id].append(t_query_gs)

    for query, label in zip(samples_metadata_queries, labels):
        for metadata in query:
            sample_id = metadata[2]['SAMPLEID']
            columns = metadata[3]
            if 'BINID' in columns:
                if not sample_id_to_g_gs[sample_id]:
                    logging.getLogger('amber').critical("Sample ID {} in {} not found in the genome binning gold standard.".format(sample_id, label))
                    exit(1)
                g_query = binning_classes.GenomeQuery(label, sample_id, options, metadata)
                g_query.gold_standard = sample_id_to_g_gs[sample_id]
                sample_id_to_g_queries_list[sample_id].append(g_query)
                options.only_taxonomic_queries = options_gs.only_taxonomic_queries = False
            if 'TAXID' in columns and not taxonomy_df.empty:
                if not sample_id_to_t_gs[sample_id]:
                    logging.getLogger('amber').critical("Sample ID {} in {} not found in the taxon binning gold standard.".format(sample_id, label))
                    exit(1)
                t_query = binning_classes.TaxonomicQuery(label, sample_id, options, metadata, taxonomy_df)
                t_query.gold_standard = sample_id_to_t_gs[sample_id]
                sample_id_to_t_queries_list[sample_id].append(t_query)
                options.only_genome_queries = options_gs.only_genome_queries = False

    return sample_id_to_g_queries_list, sample_id_to_t_queries_list, [metadata[2]['SAMPLEID'] for metadata in samples_metadata_gs]
