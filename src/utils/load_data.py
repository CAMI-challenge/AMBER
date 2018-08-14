#!/usr/bin/env python

import os
import gzip
import mimetypes
import sys
from Bio import SeqIO
import numpy as np
import warnings
import traceback
from collections import defaultdict

try:
    import exclude_genomes
    import load_ncbi_taxinfo
except ImportError:
    sys.path.append(os.path.dirname(__file__))
    try:
        import exclude_genomes
        import load_ncbi_taxinfo
    finally:
        sys.path.remove(os.path.dirname(__file__))


class NCBI:
    def __init__(self):
        self.__tax_id_to_rank = None

    @property
    def tax_id_to_rank(self):
        return self.__tax_id_to_rank

    @tax_id_to_rank.setter
    def tax_id_to_rank(self, tax_id_to_rank):
        self.__tax_id_to_rank = tax_id_to_rank


class Query:
    def __init__(self):
        self.__bins = []
        self.__bin_id_to_bin = {}
        self.__label = ""

    @property
    def bins(self):
        return self.__bins

    @property
    def label(self):
        return self.__label

    @bins.setter
    def bins(self, bins):
        self.__bins = bins

    @label.setter
    def label(self, label):
        self.__label = label

    def add_bin(self, bin):
        self.__bins.append(bin)
        self.__bin_id_to_bin[bin.id] = bin

    def get_bin_ids(self):
        return self.__bin_id_to_bin.keys()
        # return set([bin.id for bin in self.__bins])

    def get_sequence_ids(self):
        return set.union(*(bin.sequence_ids for bin in self.__bins))

    def get_bin_by_id(self, id):
        return self.__bin_id_to_bin[id]


class GenomeQuery(Query):
    def __init__(self):
        super().__init__()
        self.__sequence_id_to_bin_id = {}
        self.__map_by_completeness = False

    @property
    def sequence_id_to_bin_id(self):
        return self.__sequence_id_to_bin_id

    @property
    def map_by_completeness(self):
        return self.__map_by_completeness

    @sequence_id_to_bin_id.setter
    def sequence_id_to_bin_id(self, sequence_id_bin_id):
        (sequence_id, bin_id) = sequence_id_bin_id
        self.__sequence_id_to_bin_id[sequence_id] = bin_id

    @map_by_completeness.setter
    def map_by_completeness(self, map_by_completeness):
        self.__map_by_completeness = map_by_completeness

    def compute_true_positives(self, gold_standard):
        for bin in self.bins:
            bin.compute_true_positives(gold_standard.genome_query, self.__map_by_completeness)


class TaxonomicQuery(Query):
    def __init__(self):
        super().__init__()

    def compute_true_positives(self, gold_standard):
        for bin in self.bins:
            bin.compute_true_positives(gold_standard.taxonomic_query)


class Bin:
    sequence_id_to_length = {}

    def __init__(self, id):
        self.__id = id
        self.__sequence_ids = set()
        self.__length = 0
        self.__true_positives = 0
        self.__mapping_id = None

    @property
    def id(self):
        return self.__id

    @property
    def sequence_ids(self):
        return self.__sequence_ids

    @property
    def length(self):
        return self.__length

    @property
    def true_positives(self):
        return self.__true_positives

    @property
    def mapping_id(self):
        return self.__mapping_id

    @id.setter
    def id(self, id):
        self.__id = id

    @sequence_ids.setter
    def sequence_ids(self, sequence_ids):
        self.__sequence_ids = sequence_ids

    @length.setter
    def length(self, length):
        self.__length = length

    @true_positives.setter
    def true_positives(self, true_positives):
        self.__true_positives = true_positives

    @mapping_id.setter
    def mapping_id(self, mapping_id):
        self.__mapping_id = mapping_id

    def add_sequence_id(self, sequence_id, sequence_length):
        self.__sequence_ids.add(sequence_id)
        self.__length += sequence_length


class GenomeBin(Bin):
    def __init__(self, id):
        super().__init__(id)
        self.__genome_id_to_length = defaultdict(int)

    def compute_confusion_matrix(self, gold_standard_query):
        if self.__genome_id_to_length:
            return
        for sequence_id in self.sequence_ids:
            genome_id = gold_standard_query.sequence_id_to_bin_id[sequence_id]
            self.__genome_id_to_length[genome_id] += Bin.sequence_id_to_length[sequence_id]

    def compute_true_positives(self, gold_standard_query, map_by_completeness):
        self.compute_confusion_matrix(gold_standard_query)
        if map_by_completeness:
            max_genome_percentage = .0
            best_gs_bin = gold_standard_query.bins[0]
            for gs_bin in gold_standard_query.bins:
                genome_percentage = self.mapping_id_to_length[gs_bin.id] / gs_bin.length
                if max_genome_percentage < genome_percentage:
                    max_genome_percentage = genome_percentage
                    best_gs_bin = gs_bin
                elif max_genome_percentage == genome_percentage and gs_bin.length > best_gs_bin.length:
                    best_gs_bin = gs_bin
            self.__mapping_id = best_gs_bin.id
            self.true_positives = self.__genome_id_to_length[best_gs_bin.id]
        else:
            self.__mapping_id = max(self.__genome_id_to_length, key=self.__genome_id_to_length.get)
            self.true_positives = self.__genome_id_to_length[self.__mapping_id]


class TaxonomicBin(Bin):
    def __init__(self, id):
        super().__init__(id)
        self.__rank = None

    @property
    def rank(self):
        return self.__rank

    @rank.setter
    def rank(self, rank):
        self.__rank = rank

    @property
    def mapping_id(self):
        return self.id

    def compute_true_positives(self, gold_standard_query):
        if not self.id in gold_standard_query.get_bin_ids():
            return
        for gs_bin in gold_standard_query.bins:
            if self.id == gs_bin.id:
                commmon_seq_ids = self.sequence_ids & gs_bin.sequence_ids
                for sequence_id in commmon_seq_ids:
                    self.true_positives += Bin.sequence_id_to_length[sequence_id]


class GoldStandard:
    def __init__(self, genome_query, taxonomic_query):
        self.__genome_query = genome_query
        self.__taxonomic_query = taxonomic_query

    @property
    def genome_query(self):
        return self.__genome_query

    @property
    def taxonomic_query(self):
        return self.__taxonomic_query

    @genome_query.setter
    def genome_query(self, genome_query):
        self.__genome_query = genome_query

    @taxonomic_query.setter
    def taxonomic_query(self, taxonomic_query):
        self.__taxonomic_query = taxonomic_query


def load_unique_common(unique_common_file_path):
    genome_to_unique_common = {}
    with open(unique_common_file_path) as read_handler:
        for line in read_handler:
            genome_to_unique_common[line.split("\t")[0]] = line.split("\t")[1].strip('\n')
    return genome_to_unique_common


def load_tsv_table(stream):
    data = []
    next(stream)
    for line in stream:
        line = line.strip()
        if len(line) == 0 or line.startswith("@"):
            continue
        row_data = line.split('\t')

        mapped_genome = row_data[0]
        real_size = int(float(row_data[5]))
        predicted_size = int(float(row_data[3]))
        correctly_predicted = int(float(row_data[4]))

        if row_data[1] != "NA" and predicted_size > 0:
            precision = float(row_data[1])
        else:
            precision = np.nan
        data.append({'mapped_genome': mapped_genome, 'precision': precision, 'recall': float(row_data[2]),
                     'predicted_size': predicted_size, 'correctly_predicted': correctly_predicted, 'real_size': real_size})
    return data


def read_lengths_from_fastx_file(fastx_file):
    """

    @param fastx_file: file path
    @type fastx_file: str
    @rtype: dict[str, int]
    """
    file_type = mimetypes.guess_type(fastx_file)[1]
    if file_type == 'gzip':
        f = gzip.open(fastx_file, "rt")
    elif not file_type:
        f = open(fastx_file, "rt")
    else:
        raise RuntimeError("Unknown type of file: '{}".format(fastx_file))

    length = {}
    if os.path.getsize(fastx_file) == 0:
        return length

    file_format = None
    line = f.readline()
    if line.startswith('@'):
        file_format = "fastq"
    elif line.startswith(">"):
        file_format = "fasta"
    f.seek(0)

    if not file_format:
        raise RuntimeError("Invalid sequence file: '{}".format(fastx_file))

    for seq_record in SeqIO.parse(f, file_format):
        length[seq_record.id] = len(seq_record.seq)

    f.close()

    return length


def get_genome_mapping_without_lenghts(mapping_file, remove_genomes_file=None, keyword=None):
    gold_standard = Query()
    gold_standard.bin_id_to_list_of_sequence_ids = {}
    gold_standard.sequence_id_to_genome_id = {}

    filtering_genomes_to_keyword = {}
    if remove_genomes_file:
        filtering_genomes_to_keyword = load_unique_common(remove_genomes_file)

    with open(mapping_file, 'r') as read_handler:
        try:
            for sequence_id, genome_id, length in read_binning_file(read_handler):
                if genome_id in filtering_genomes_to_keyword and (not keyword or filtering_genomes_to_keyword[genome_id] == keyword):
                    continue
                gold_standard.sequence_id_to_genome_id[sequence_id] = genome_id
                if genome_id not in gold_standard.bin_id_to_list_of_sequence_ids:
                    gold_standard.bin_id_to_list_of_sequence_ids[genome_id] = []
                gold_standard.bin_id_to_list_of_sequence_ids[genome_id].append(sequence_id)
        except:
            exit("Error. File {} is malformed.".format(mapping_file))

    if len(gold_standard.bin_id_to_list_of_sequence_ids) == 0:
        exit('All bins of the gold standard have been filtered out due to option --remove_genomes.')

    return gold_standard


def read_header(input_stream):
    """
    @Version:0.9.1
    @SampleID:RH_S1

    @@SEQUENCEID    BINID   TAXID
    """
    header = {}
    column_names = {}
    for line in input_stream:
        if len(line.strip()) == 0 or line.startswith("#"):
            continue
        line = line.rstrip('\n')
        if line.startswith("@@"):
            for index, column_name in enumerate(line[2:].split('\t')):
                column_names[column_name] = index
            return header, column_names
        if line.startswith("@"):
            key, value = line[1:].split(':', 1)
            header[key] = value.strip()


def read_rows(input_stream, index_seq_id, index_bin_id, index_tax_id, index_length, is_gs):
    for line in input_stream:
        if len(line.strip()) == 0 or line.startswith("#"):
            continue
        line = line.rstrip('\n')
        row_data = line.split('\t')

        try:
            seq_id = row_data[index_seq_id]
        except:
            print("Value in column SEQUENCEID could not be read.")
            raise

        if index_bin_id:
            try:
                bin_id = row_data[index_bin_id]
            except:
                print("Value in column BINID could not be read.")
                raise
        else:
            bin_id = None

        if index_tax_id:
            try:
                tax_id = row_data[index_tax_id]
            except:
                print("Value in column TAXID could not be read.")
                raise
        else:
            tax_id = None

        if is_gs and index_length is not None:
            try:
                length = int(row_data[index_length])
            except:
                print("Value in column _LENGTH could not be read. Please provide a value or remove column altogether (and provide a FASTA or FASTQ file instead - see README).")
                raise
            yield seq_id, bin_id, tax_id, length
        else:
            yield seq_id, bin_id, tax_id, int(0)


def is_length_column_available(input_stream):
    header, column_names = read_header(input_stream)
    input_stream.seek(0)
    return "_LENGTH" in column_names


def get_column_indices(input_stream):
    """

    :param input_stream: 
    :return: 
    """
    header, column_names = read_header(input_stream)
    if "SEQUENCEID" not in column_names:
        raise RuntimeError("Column not found: {}".format("SEQUENCEID"))
    if "BINID" not in column_names and "TAXID" not in column_names:
        raise RuntimeError("Column not found: {}".format("BINID/TAXID"))
    index_seq_id = column_names["SEQUENCEID"]

    index_bin_id = None
    index_tax_id = None
    if "BINID" in column_names:
        index_bin_id = column_names["BINID"]
    if "TAXID" in column_names:
        index_tax_id = column_names["TAXID"]

    index_length = None
    if "_LENGTH" in column_names:
        index_length = column_names["_LENGTH"]
    return index_seq_id, index_bin_id, index_tax_id, index_length


def read_binning_file(input_stream, is_gs):
    index_seq_id, index_bin_id, index_tax_id, index_length = get_column_indices(input_stream)
    return read_rows(input_stream, index_seq_id, index_bin_id, index_tax_id, index_length, is_gs)


def open_query(file_path_query, is_gs, fastx_file, tax_id_to_parent, tax_id_to_rank, gold_standard, min_length=0):
    g_query = GenomeQuery()
    t_query = TaxonomicQuery()
    is_length_column_av = False

    with open(file_path_query) as read_handler:
        if is_gs and not Bin.sequence_id_to_length:
            is_length_column_av = is_length_column_available(read_handler)
            if not is_length_column_av:
                if not fastx_file:
                    exit("Sequences length could not be determined. Please provide a FASTA or FASTQ file using option -f or add column _LENGTH to gold standard.")
                Bin.sequence_id_to_length = read_lengths_from_fastx_file(fastx_file)

        if gold_standard and gold_standard.genome_query:
            gs_sequence_ids = gold_standard.genome_query.get_sequence_ids()

        try:
            for sequence_id, bin_id, tax_id, length in read_binning_file(read_handler, is_gs):
                if is_gs and is_length_column_av:
                    Bin.sequence_id_to_length[sequence_id] = length

                if bin_id:
                    if gold_standard and sequence_id not in gs_sequence_ids:
                        continue
                    if Bin.sequence_id_to_length[sequence_id] < min_length:
                        continue
                    if bin_id not in g_query.get_bin_ids():
                        bin = GenomeBin(bin_id)
                        g_query.add_bin(bin)
                    else:
                        bin = g_query.get_bin_by_id(bin_id)
                    g_query.sequence_id_to_bin_id = (sequence_id, bin_id)
                    bin.add_sequence_id(sequence_id, Bin.sequence_id_to_length[sequence_id])

                if tax_id:
                    # TODO: check for every rank
                    # if t_gold_standard and sequence_id not in t_gold_standard.sequence_id_to_bin_id:
                    #     continue
                    if not tax_id_to_parent:
                        warnings.warn("Taxonomic binning cannot be assessed. Please provide an NCBI nodes file using option --ncbi_nodes_file.", Warning)
                    tax_id_path = load_ncbi_taxinfo.get_id_path(tax_id, tax_id_to_parent, tax_id_to_rank)
                    for tax_id in tax_id_path:
                        if not tax_id: # tax_id may be empty
                            continue
                        if tax_id not in t_query.get_bin_ids():
                            bin = TaxonomicBin(tax_id)
                            t_query.add_bin(bin)
                        else:
                            bin = t_query.get_bin_by_id(tax_id)

                        bin.add_sequence_id(sequence_id, Bin.sequence_id_to_length[sequence_id])

                        # rank = tax_id_to_rank[tax_id]
                        # if rank not in t_query.rank_to_sequence_id_to_bin_id:
                        #     t_query.rank_to_sequence_id_to_bin_id[rank] = {}
                        # t_query.rank_to_sequence_id_to_bin_id[rank][sequence_id] = tax_id
        except BaseException as e:
            traceback.print_exc()
            exit("Error. File {} is malformed. {}".format(file_path_query, e))

    if not g_query.bins:
        g_query = None
    if not t_query.bins:
        t_query = None
    return g_query, t_query

