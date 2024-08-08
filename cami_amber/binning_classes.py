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
import logging
import itertools
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')
import os
import math
from abc import ABC, abstractmethod
from collections import defaultdict
from collections import OrderedDict
from cami_amber.utils import labels as utils_labels
from cami_amber.utils import load_data
from cami_amber.utils import load_ncbi_taxinfo
from cami_amber.utils import ProfilingTools as pf
from cami_amber import unifrac_distance as uf
from cami_amber import precision_recall_per_bin
from cami_amber import plots


class Metrics:
    def __init__(self):
        self.__percentage_of_assigned_bps = .0
        self.__percentage_of_assigned_seqs = .0
        self.__accuracy_bp = .0
        self.__accuracy_seq = .0
        self.__rand_index_bp = .0
        self.__adjusted_rand_index_bp = .0
        self.__rand_index_seq = .0
        self.__adjusted_rand_index_seq = .0
        self.__precision_avg_bp = .0
        self.__precision_avg_bp_var = .0
        self.__precision_avg_bp_sem = .0
        self.__precision_avg_seq = .0
        self.__precision_avg_seq_sem = .0
        self.__precision_weighted_bp = .0
        self.__precision_weighted_seq = .0
        self.__recall_avg_bp = .0
        self.__recall_avg_bp_cami1 = .0
        self.__recall_avg_bp_var = .0
        self.__recall_avg_bp_var_cami1 = .0
        self.__recall_avg_bp_sem = .0
        self.__recall_avg_bp_sem_cami1 = .0
        self.__recall_avg_seq = .0
        self.__recall_avg_seq_cami1 = .0
        self.__recall_avg_seq_sem = .0
        self.__recall_avg_seq_sem_cami1 = .0
        self.__recall_weighted_bp = .0
        self.__recall_weighted_seq = .0
        self.__f1_score_bp = .0
        self.__f1_score_seq = .0
        self.__f1_score_bp_cami1 = .0
        self.__f1_score_seq_cami1 = .0
        self.__unifrac_bp = None
        self.__unifrac_seq = None

    @property
    def percentage_of_assigned_bps(self):
        return self.__percentage_of_assigned_bps

    @property
    def percentage_of_assigned_seqs(self):
        return self.__percentage_of_assigned_seqs

    @property
    def accuracy_bp(self):
        return self.__accuracy_bp

    @property
    def accuracy_seq(self):
        return self.__accuracy_seq

    @property
    def rand_index_bp(self):
        return self.__rand_index_bp

    @property
    def adjusted_rand_index_bp(self):
        return self.__adjusted_rand_index_bp

    @property
    def rand_index_seq(self):
        return self.__rand_index_seq

    @property
    def adjusted_rand_index_seq(self):
        return self.__adjusted_rand_index_seq

    @property
    def precision_avg_bp(self):
        return self.__precision_avg_bp

    @property
    def precision_avg_bp_var(self):
        return self.__precision_avg_bp_var

    @property
    def precision_avg_bp_sem(self):
        return self.__precision_avg_bp_sem

    @property
    def precision_avg_seq(self):
        return self.__precision_avg_seq

    @property
    def precision_avg_seq_sem(self):
        return self.__precision_avg_seq_sem

    @property
    def precision_weighted_bp(self):
        return self.__precision_weighted_bp

    @property
    def precision_weighted_seq(self):
        return self.__precision_weighted_seq

    @property
    def recall_avg_bp(self):
        return self.__recall_avg_bp

    @property
    def recall_avg_bp_cami1(self):
        return self.__recall_avg_bp_cami1

    @property
    def recall_avg_bp_var(self):
        return self.__recall_avg_bp_var

    @property
    def recall_avg_bp_var_cami1(self):
        return self.__recall_avg_bp_var_cami1

    @property
    def recall_avg_bp_sem(self):
        return self.__recall_avg_bp_sem

    @property
    def recall_avg_bp_sem_cami1(self):
        return self.__recall_avg_bp_sem_cami1

    @property
    def recall_avg_seq(self):
        return self.__recall_avg_seq

    @property
    def recall_avg_seq_cami1(self):
        return self.__recall_avg_seq_cami1

    @property
    def recall_avg_seq_sem(self):
        return self.__recall_avg_seq_sem

    @property
    def recall_avg_seq_sem_cami1(self):
        return self.__recall_avg_seq_sem_cami1

    @property
    def recall_weighted_bp(self):
        return self.__recall_weighted_bp

    @property
    def recall_weighted_seq(self):
        return self.__recall_weighted_seq

    @property
    def unifrac_bp(self):
        return self.__unifrac_bp

    @property
    def unifrac_seq(self):
        return self.__unifrac_seq

    @percentage_of_assigned_bps.setter
    def percentage_of_assigned_bps(self, percentage_of_assigned_bps):
        self.__percentage_of_assigned_bps = percentage_of_assigned_bps

    @percentage_of_assigned_seqs.setter
    def percentage_of_assigned_seqs(self, percentage_of_assigned_seqs):
        self.__percentage_of_assigned_seqs = percentage_of_assigned_seqs

    @accuracy_bp.setter
    def accuracy_bp(self, accuracy_bp):
        self.__accuracy_bp = accuracy_bp

    @accuracy_seq.setter
    def accuracy_seq(self, accuracy_seq):
        self.__accuracy_seq = accuracy_seq

    @rand_index_bp.setter
    def rand_index_bp(self, rand_index_bp):
        self.__rand_index_bp = rand_index_bp

    @adjusted_rand_index_bp.setter
    def adjusted_rand_index_bp(self, adjusted_rand_index_bp):
        self.__adjusted_rand_index_bp = adjusted_rand_index_bp

    @rand_index_seq.setter
    def rand_index_seq(self, rand_index_seq):
        self.__rand_index_seq = rand_index_seq

    @adjusted_rand_index_seq.setter
    def adjusted_rand_index_seq(self, adjusted_rand_index_seq):
        self.__adjusted_rand_index_seq = adjusted_rand_index_seq

    @precision_avg_bp.setter
    def precision_avg_bp(self, precision_avg_bp):
        self.__precision_avg_bp = precision_avg_bp

    @precision_avg_bp_var.setter
    def precision_avg_bp_var(self, precision_avg_bp_var):
        self.__precision_avg_bp_var = precision_avg_bp_var

    @precision_avg_bp_sem.setter
    def precision_avg_bp_sem(self, precision_avg_bp_sem):
        self.__precision_avg_bp_sem = precision_avg_bp_sem

    @precision_avg_seq.setter
    def precision_avg_seq(self, precision_avg_seq):
        self.__precision_avg_seq = precision_avg_seq

    @precision_avg_seq_sem.setter
    def precision_avg_seq_sem(self, precision_avg_seq_sem):
        self.__precision_avg_seq_sem = precision_avg_seq_sem

    @precision_weighted_bp.setter
    def precision_weighted_bp(self, precision_weighted_bp):
        self.__precision_weighted_bp = precision_weighted_bp

    @precision_weighted_seq.setter
    def precision_weighted_seq(self, precision_weighted_seq):
        self.__precision_weighted_seq = precision_weighted_seq

    @recall_avg_bp.setter
    def recall_avg_bp(self, recall_avg_bp):
        self.__recall_avg_bp = recall_avg_bp

    @recall_avg_bp_cami1.setter
    def recall_avg_bp_cami1(self, recall_avg_bp_cami1):
        self.__recall_avg_bp_cami1 = recall_avg_bp_cami1

    @recall_avg_bp_var.setter
    def recall_avg_bp_var(self, recall_avg_bp_var):
        self.__recall_avg_bp_var = recall_avg_bp_var

    @recall_avg_bp_var_cami1.setter
    def recall_avg_bp_var_cami1(self, recall_avg_bp_var_cami1):
        self.__recall_avg_bp_var_cami1 = recall_avg_bp_var_cami1

    @recall_avg_bp_sem.setter
    def recall_avg_bp_sem(self, recall_avg_bp_sem):
        self.__recall_avg_bp_sem = recall_avg_bp_sem

    @recall_avg_bp_sem_cami1.setter
    def recall_avg_bp_sem_cami1(self, recall_avg_bp_sem_cami1):
        self.__recall_avg_bp_sem_cami1 = recall_avg_bp_sem_cami1

    @recall_avg_seq.setter
    def recall_avg_seq(self, recall_avg_seq):
        self.__recall_avg_seq = recall_avg_seq

    @recall_avg_seq_cami1.setter
    def recall_avg_seq_cami1(self, recall_avg_seq_cami1):
        self.__recall_avg_seq_cami1 = recall_avg_seq_cami1

    @recall_avg_seq_sem.setter
    def recall_avg_seq_sem(self, recall_avg_seq_sem):
        self.__recall_avg_seq_sem = recall_avg_seq_sem

    @recall_avg_seq_sem_cami1.setter
    def recall_avg_seq_sem_cami1(self, recall_avg_seq_sem_cami1):
        self.__recall_avg_seq_sem_cami1 = recall_avg_seq_sem_cami1

    @recall_weighted_bp.setter
    def recall_weighted_bp(self, recall_weighted_bp):
        self.__recall_weighted_bp = recall_weighted_bp

    @recall_weighted_seq.setter
    def recall_weighted_seq(self, recall_weighted_seq):
        self.__recall_weighted_seq = recall_weighted_seq

    @unifrac_bp.setter
    def unifrac_bp(self, unifrac_bp):
        self.__unifrac_bp = unifrac_bp

    @unifrac_seq.setter
    def unifrac_seq(self, unifrac_seq):
        self.__unifrac_seq = unifrac_seq

    @staticmethod
    def compute_rand_index(confusion_df, col_name, gs_col_name, field):
        def choose2(n):
            return math.comb(n, 2)

        bin_mapping_comb = sum(confusion_df[field].apply(choose2))
        bin_comb = sum(confusion_df.groupby(col_name).agg({field: 'sum'})[field].apply(choose2))
        mapping_comb = sum(confusion_df.groupby(gs_col_name).agg({field: 'sum'})[field].apply(choose2))
        num_bp_comb = choose2(sum(confusion_df[field]))

        rand_index = ((num_bp_comb - bin_comb - mapping_comb + 2 * bin_mapping_comb) / num_bp_comb) if num_bp_comb != 0 else .0

        temp = (bin_comb * (mapping_comb / num_bp_comb)) if num_bp_comb != 0 else .0
        ret = bin_mapping_comb - temp
        denominator = (((bin_comb + mapping_comb) / 2.0) - temp)
        adjusted_rand_index = (ret / denominator) if denominator != 0 else .0

        return rand_index, adjusted_rand_index

    def get_ordered_dict(self):
        def f1_score(metric1, metric2):
            if metric1 + metric2 > 0:
                return 2 * metric1 * metric2 / (metric1 + metric2)
            else:
                return np.nan

        metrics_dict = {metric.split("__")[1]: value for metric, value in self.__dict__.items()}
        metrics_dict['f1_score_bp'] = f1_score(self.__precision_avg_bp, self.__recall_avg_bp)
        metrics_dict['f1_score_seq'] = f1_score(self.__precision_avg_seq, self.__recall_avg_seq)
        metrics_dict['f1_score_bp_cami1'] = f1_score(self.__precision_avg_bp, self.__recall_avg_bp_cami1)
        metrics_dict['f1_score_seq_cami1'] = f1_score(self.__precision_avg_seq, self.__recall_avg_seq_cami1)
        metrics_dict['f1_score_per_bp'] = f1_score(self.__precision_weighted_bp, self.__recall_weighted_bp)
        metrics_dict['f1_score_per_seq'] = f1_score(self.__precision_weighted_seq, self.__recall_weighted_seq)
        metrics_dict['misclassification_bp'] = 1 - self.__precision_weighted_bp
        metrics_dict['misclassification_seq'] = 1 - self.__precision_weighted_seq

        return OrderedDict([(utils_labels.TOOL, None),
                            (utils_labels.BINNING_TYPE, None),
                            (utils_labels.SAMPLE, None),
                            (utils_labels.RANK, None)] +
                            [(k, [v]) for k, v in metrics_dict.items()])


class Query(ABC):
    def __init__(self, label, sample_id, options, metadata, is_gs=False):
        self.__label = label
        self.__sample_id = sample_id
        self.__gold_standard = None
        self.__precision_df = pd.DataFrame()
        self.__recall_df = pd.DataFrame()
        self.__heatmap_sdf = None
        self.__metrics = None
        self.__metrics_filtered = None
        self.__options = options
        self.__metadata = metadata
        self.__is_gs = is_gs

    @property
    def label(self):
        return self.__label

    @property
    def sample_id(self):
        return self.__sample_id

    @property
    def gold_standard(self):
        return self.__gold_standard

    @property
    def precision_df(self):
        return self.__precision_df

    @property
    def recall_df(self):
        return self.__recall_df

    @property
    def heatmap_sdf(self):
        return self.__heatmap_sdf

    @property
    def metrics(self):
        return self.__metrics

    @property
    def metrics_filtered(self):
        return self.__metrics_filtered

    @property
    def options(self):
        return self.__options

    @property
    def metadata(self):
        return self.__metadata

    @property
    def is_gs(self):
        return self.__is_gs

    @label.setter
    def label(self, label):
        self.__label = label

    @sample_id.setter
    def sample_id(self, sample_id):
        self.__sample_id = sample_id

    @gold_standard.setter
    def gold_standard(self, gold_standard):
        self.__gold_standard = gold_standard

    @precision_df.setter
    def precision_df(self, precision_df):
        self.__precision_df = precision_df

    @recall_df.setter
    def recall_df(self, recall_df):
        self.__recall_df = recall_df

    @heatmap_sdf.setter
    def heatmap_sdf(self, heatmap_sdf):
        self.__heatmap_sdf = heatmap_sdf

    @metrics.setter
    def metrics(self, metrics):
        self.__metrics = metrics

    @metrics_filtered.setter
    def metrics_filtered(self, metrics_filtered):
        self.__metrics_filtered = metrics_filtered

    @options.setter
    def options(self, options):
        self.__options = options

    @metadata.setter
    def metadata(self, metadata):
        self.__metadata = metadata

    @abstractmethod
    def compute_metrics(self):
        pass

    def compute_unifrac(self, all_bins):
        return None, None

    def plot(self):
        return


class GenomeQuery(Query):
    binning_type = 'genome'

    def __init__(self, label, sample_id, options, metadata, is_gs=False):
        super().__init__(label, sample_id, options, metadata, is_gs)
        self.__recall_df_cami1 = pd.DataFrame()
        self.metrics = Metrics()

    @property
    def df(self):
        query_df = load_data.load_sample(self.metadata).drop_duplicates(['SEQUENCEID', 'BINID'])
        if self.options.min_length:
            query_df = query_df[query_df['LENGTH'] >= self.options.min_length]
        if self.is_gs:
            return query_df[['SEQUENCEID', 'BINID', 'LENGTH']]
        else:
            return query_df[['SEQUENCEID', 'BINID']]

    @property
    def recall_df_cami1(self):
        return self.__recall_df_cami1

    @recall_df_cami1.setter
    def recall_df_cami1(self, recall_df_cami1):
        self.__recall_df_cami1 = recall_df_cami1

    def get_metrics_df(self):
        metrics_dict = self.metrics.get_ordered_dict()
        metrics_dict[utils_labels.TOOL] = self.label
        metrics_dict[utils_labels.BINNING_TYPE] = self.binning_type
        metrics_dict[utils_labels.RANK] = 'NA'
        return pd.DataFrame(metrics_dict)

    def compute_metrics(self, gs_df):
        if self.label == utils_labels.GS and (self.options.only_taxonomic_queries or self.options.skip_gs):
            return False
        logging.getLogger('amber').info('Evaluating {}, sample {}, genome binning'.format(self.label, self.sample_id))

        query_df = self.df
        condition = query_df['SEQUENCEID'].isin(gs_df['SEQUENCEID'])
        if ~condition.all():
            logging.getLogger('amber').warning("{} sequences in {} not found in the gold standard.".format(query_df[~condition]['SEQUENCEID'].nunique(), self.label))
            query_df = query_df[condition]

        query_w_length = pd.merge(query_df, gs_df.drop_duplicates('SEQUENCEID'), on='SEQUENCEID', sort=False)

        query_w_length_no_dups = query_w_length.drop_duplicates('SEQUENCEID')
        gs_df_no_dups = gs_df.drop_duplicates('SEQUENCEID')
        self.metrics.percentage_of_assigned_bps = query_w_length_no_dups['seq_length'].sum() / gs_df_no_dups['seq_length'].sum()
        self.metrics.percentage_of_assigned_seqs = query_w_length_no_dups.shape[0] / gs_df_no_dups['SEQUENCEID'].shape[0]
        del query_w_length_no_dups

        # confusion table possibly with the same sequences in multiple bins
        query_w_length_mult_seqs = query_df.reset_index().merge(gs_df, on='SEQUENCEID', sort=False)

        if query_w_length.shape[0] < query_w_length_mult_seqs.shape[0]:
            query_w_length_mult_seqs.drop_duplicates(['index', 'genome_id'], inplace=True)
            confusion_df = query_w_length_mult_seqs.groupby(['BINID', 'genome_id'], sort=False).agg({'seq_length': 'sum', 'SEQUENCEID': 'count'}).rename(columns={'seq_length': 'genome_length', 'SEQUENCEID': 'genome_seq_counts'})

            most_abundant_genome_df = confusion_df.loc[confusion_df.groupby('BINID', sort=False)['genome_length'].idxmax()]
            most_abundant_genome_df = most_abundant_genome_df.reset_index()[['BINID', 'genome_id']]

            matching_genomes_df = pd.merge(query_w_length_mult_seqs, most_abundant_genome_df, on=['BINID', 'genome_id']).set_index('index')
            query_w_length_mult_seqs.set_index('index', inplace=True)
            difference_df = query_w_length_mult_seqs.drop(matching_genomes_df.index).groupby(['index'], sort=False).first()
            query_w_length = pd.concat([matching_genomes_df, difference_df])

            # Modify gs such that multiple binnings of the same sequence are not required
            matching_genomes_df = pd.merge(gs_df, query_w_length[['SEQUENCEID', 'genome_id']], on=['SEQUENCEID', 'genome_id'])
            matching_genomes_df = matching_genomes_df[['SEQUENCEID', 'genome_id', 'seq_length']].drop_duplicates(['SEQUENCEID', 'genome_id'])
            condition = gs_df_no_dups['SEQUENCEID'].isin(matching_genomes_df['SEQUENCEID'])
            difference_df = gs_df_no_dups[~condition]
            gs_df = pd.concat([difference_df, matching_genomes_df])

            # query_w_length_mult_seqs.reset_index(inplace=True)
            # query_w_length_mult_seqs = pd.merge(query_w_length_mult_seqs, most_abundant_genome_df, on=['BINID'])
            # grouped = query_w_length_mult_seqs.groupby(['index'], sort=False, as_index=False)
            # query_w_length = grouped.apply(lambda x: x[x['genome_id_x'] == x['genome_id_y'] if any(x['genome_id_x'] == x['genome_id_y']) else len(x) * [True]])
            # query_w_length = query_w_length.groupby(['index'], sort=False).first().drop(columns='genome_id_y').rename(columns={'genome_id_x': 'genome_id'})

        del query_w_length_mult_seqs

        confusion_df = query_w_length.groupby(['BINID', 'genome_id'], sort=False).agg({'seq_length': 'sum', 'SEQUENCEID': 'count'}).rename(columns={'seq_length': 'genome_length', 'SEQUENCEID': 'genome_seq_counts'})

        self.metrics.rand_index_bp, self.metrics.adjusted_rand_index_bp = Metrics.compute_rand_index(confusion_df, 'BINID', 'genome_id', 'genome_length')
        self.metrics.rand_index_seq, self.metrics.adjusted_rand_index_seq = Metrics.compute_rand_index(confusion_df, 'BINID', 'genome_id', 'genome_seq_counts')

        most_abundant_genome_df = confusion_df.loc[confusion_df.groupby('BINID', sort=False)['genome_length'].idxmax()].reset_index().set_index('BINID')

        query_w_length['seq_length_mean'] = query_w_length['seq_length']

        precision_df = query_w_length.groupby('BINID', sort=False).agg({'seq_length': 'sum', 'seq_length_mean': 'mean', 'SEQUENCEID': 'count'}).rename(columns={'seq_length': 'total_length', 'SEQUENCEID': 'total_seq_counts'})
        precision_df = pd.merge(precision_df, most_abundant_genome_df, on='BINID')
        precision_df.rename(columns={'genome_length': 'tp_length', 'genome_seq_counts': 'tp_seq_counts'}, inplace=True)
        precision_df['precision_bp'] = precision_df['tp_length'] / precision_df['total_length']
        precision_df['precision_seq'] = precision_df['tp_seq_counts'] / precision_df['total_seq_counts']

        if self.options.filter_tail_percentage:
            precision_df['total_length_pct'] = precision_df['total_length'] / precision_df['total_length'].sum()
            precision_df.sort_values(by='total_length', inplace=True)
            precision_df['cumsum_length_pct'] = precision_df['total_length_pct'].cumsum(axis=0)
            precision_df['precision_bp'].mask(precision_df['cumsum_length_pct'] <= self.options.filter_tail_percentage / 100, inplace=True)
            precision_df['precision_seq'].mask(precision_df['precision_bp'].isna(), inplace=True)
            precision_df.drop(columns=['cumsum_length_pct', 'total_length_pct'], inplace=True)
        if self.options.genome_to_unique_common:
            precision_df = precision_df[~precision_df['genome_id'].isin(self.options.genome_to_unique_common)]

        self.metrics.precision_avg_bp = precision_df['precision_bp'].mean()
        self.metrics.precision_avg_bp_sem = precision_df['precision_bp'].sem()
        self.metrics.precision_avg_bp_var = precision_df['precision_bp'].var()
        self.metrics.precision_avg_seq = precision_df['precision_seq'].mean()
        self.metrics.precision_avg_seq_sem = precision_df['precision_seq'].sem()

        def safe_divide(x, y):
            return x / y if y else np.nan

        self.metrics.precision_weighted_bp = safe_divide(precision_df['tp_length'].sum(), precision_df['total_length'].sum())
        self.metrics.precision_weighted_seq = safe_divide(precision_df['tp_seq_counts'].sum(), precision_df['total_seq_counts'].sum())

        genome_sizes_df = gs_df.groupby('genome_id', sort=False).agg({'seq_length': 'sum', 'SEQUENCEID': 'count'}).rename(columns={'seq_length': 'length_gs', 'SEQUENCEID': 'seq_counts_gs'})
        precision_df = precision_df.reset_index().join(genome_sizes_df, on='genome_id', how='left', sort=False).set_index('BINID')
        precision_df['recall_bp'] = precision_df['tp_length'] / precision_df['length_gs']
        precision_df['recall_seq'] = precision_df['tp_seq_counts'] / precision_df['seq_counts_gs']
        precision_df['rank'] = 'NA'

        recall_df = confusion_df.loc[confusion_df.groupby('genome_id', sort=False)['genome_length'].idxmax()]
        recall_df = recall_df.reset_index().join(genome_sizes_df, on='genome_id', how='right', sort=False).set_index('BINID')
        recall_df.fillna({'genome_length': 0, 'genome_seq_counts': 0}, inplace=True)
        recall_df['recall_bp'] = recall_df['genome_length'] / recall_df['length_gs']
        recall_df['recall_seq'] = recall_df['genome_seq_counts'] / recall_df['seq_counts_gs']

        recall_df = recall_df.join(precision_df[['total_length', 'seq_length_mean']], how='left', sort=False)

        if self.options.genome_to_unique_common:
            recall_df = recall_df[~recall_df['genome_id'].isin(self.options.genome_to_unique_common)]

        self.metrics.recall_avg_bp = recall_df['recall_bp'].mean()
        self.metrics.recall_avg_bp_var = recall_df['recall_bp'].var()
        self.metrics.recall_avg_bp_sem = recall_df['recall_bp'].sem()
        self.metrics.recall_avg_seq = recall_df['recall_seq'].mean()
        self.metrics.recall_avg_seq_sem = recall_df['recall_seq'].sem()
        self.metrics.recall_weighted_bp = recall_df['genome_length'].sum() / recall_df['length_gs'].sum()
        self.metrics.recall_weighted_seq = recall_df['genome_seq_counts'].sum() / recall_df['seq_counts_gs'].sum()

        # Compute recall as in CAMI 1
        unmapped_genomes = set(gs_df['genome_id'].unique()) - set(precision_df['genome_id'].unique())
        if self.options.genome_to_unique_common:
            unmapped_genomes -= set(self.options.genome_to_unique_common)
        num_unmapped_genomes = len(unmapped_genomes)
        prec_copy = precision_df.reset_index()
        if num_unmapped_genomes:
            numeric_cols = prec_copy.select_dtypes('number').columns
            prec_copy = prec_copy.reindex(prec_copy.index.tolist() + list(range(len(prec_copy), len(prec_copy) + num_unmapped_genomes))).fillna(dict.fromkeys(numeric_cols, 0))
        self.metrics.recall_avg_bp_cami1 = prec_copy['recall_bp'].mean()
        self.metrics.recall_avg_seq_cami1 = prec_copy['recall_seq'].mean()
        self.metrics.recall_avg_bp_sem_cami1 = prec_copy['recall_bp'].sem()
        self.metrics.recall_avg_seq_sem_cami1 = prec_copy['recall_seq'].sem()
        self.metrics.recall_avg_bp_var_cami1 = prec_copy['recall_bp'].var()
        self.recall_df_cami1 = prec_copy
        # End Compute recall as in CAMI 1

        self.metrics.accuracy_bp = precision_df['tp_length'].sum() / recall_df['length_gs'].sum()
        self.metrics.accuracy_seq = precision_df['tp_seq_counts'].sum() / recall_df['seq_counts_gs'].sum()

        self.precision_df = precision_df.sort_values(by=['recall_bp'], axis=0, ascending=False)
        self.recall_df = recall_df

        self.heatmap_sdf = precision_recall_per_bin.transform_confusion_matrix2(query_w_length, confusion_df, precision_df, gs_df, log_scale=True)

        return True

    @staticmethod
    def calc_num_recovered_genomes(pd_bins, min_completeness, max_contamination):
        counts_list = []
        for (sample_id, tool), pd_group in pd_bins.groupby(['sample_id', utils_labels.TOOL]):
            for x in itertools.product(min_completeness, max_contamination):
                count = pd_group[(pd_group['recall_bp'] > x[0]) & (pd_group['precision_bp'] > (1 - x[1]))].shape[0]
                counts_list.append((sample_id, tool, '> ' + str(int(x[0] * 100)) + '% completeness', '< ' + str(int(x[1] * 100)) + '%', count))

        pd_counts = pd.DataFrame(counts_list, columns=[utils_labels.SAMPLE, utils_labels.TOOL, 'Completeness', 'Contamination', 'count'])
        pd_counts = pd.pivot_table(pd_counts, values='count', index=[utils_labels.SAMPLE, utils_labels.TOOL, 'Contamination'], columns=['Completeness']).reset_index()
        return pd_counts

    def plot_precision_vs_bin_size(self):
        fig, axs = plt.subplots(figsize=(5, 4.5))
        df_sorted = self.precision_df[['total_length', 'precision_bp']].sort_values(by=['total_length'])

        axs.scatter(np.log(df_sorted['total_length']), df_sorted['precision_bp'], marker='o')
        window = int(df_sorted.shape[0] / 50) if df_sorted.shape[0] > 100 else int(df_sorted.shape[0] / 10)
        rolling_mean = df_sorted['precision_bp'].rolling(window=window, min_periods=int(window/2)).mean()
        axs.plot(np.log(df_sorted['total_length']), rolling_mean, color='orange')

        axs.set_xlim([None, np.log(df_sorted['total_length'].max())])
        axs.set_ylim([0.0, 1.0])
        axs.set_title(self.label, fontsize=12)
        plt.ylabel('Purity per bin (%)', fontsize=12)
        plt.xlabel('Bin size [log(# bp)]', fontsize=12)
        fig.savefig(os.path.join(self.options.output_dir, 'genome', self.label, 'purity_vs_bin_size_' + self.sample_id + '.png'), dpi=200, format='png', bbox_inches='tight')
        fig.savefig(os.path.join(self.options.output_dir, 'genome', self.label, 'purity_vs_bin_size_' + self.sample_id + '.pdf'), dpi=200, format='pdf', bbox_inches='tight')
        plt.close(fig)

    def plot_recall_vs_genome_size(self):
        fig, axs = plt.subplots(figsize=(5, 4.5))
        df_sorted = self.recall_df[['total_length', 'recall_bp']].sort_values(by=['total_length'])

        axs.scatter(np.log(df_sorted['total_length']), df_sorted['recall_bp'], marker='o')
        window = int(df_sorted.shape[0] / 50) if df_sorted.shape[0] > 100 else int(df_sorted.shape[0] / 10)
        rolling_mean = df_sorted['recall_bp'].rolling(window=window, min_periods=int(window / 2)).mean()
        axs.plot(np.log(df_sorted['total_length']), rolling_mean, color='orange')

        axs.set_xlim([None, np.log(self.recall_df['total_length'].max())])
        axs.set_ylim([0.0, 1.0])
        axs.set_title(self.label, fontsize=12)
        plt.ylabel('Completeness per genome (%)', fontsize=12)
        plt.xlabel('Genome size [log(# bp)]', fontsize=12)
        fig.savefig(os.path.join(self.options.output_dir, 'genome', self.label, 'completeness_vs_genome_size_' + self.sample_id + '.png'), dpi=200, format='png', bbox_inches='tight')
        fig.savefig(os.path.join(self.options.output_dir, 'genome', self.label, 'completeness_vs_genome_size_' + self.sample_id + '.pdf'), dpi=200, format='pdf', bbox_inches='tight')
        plt.close(fig)

    def plot_heat_maps(self):
        if self.label == utils_labels.GS:
            return
        plots.plot_heatmap(self.heatmap_sdf, self.sample_id, self.options.output_dir, self.label, log_scale=True)

    def plot(self):
        if self.label == utils_labels.GS:
            return
        self.plot_precision_vs_bin_size()
        self.plot_recall_vs_genome_size()
        self.plot_heat_maps()


class TaxonomicQuery(Query):
    binning_type = 'taxonomic'

    def __init__(self, label, sample_id, options, metadata, taxonomy_df, is_gs=False):
        super().__init__(label, sample_id, options, metadata, is_gs)
        self.metrics = defaultdict()
        if self.options.filter_tail_percentage:
            self.metrics_filtered = defaultdict()
        self.__profile = None
        self.__profile_filtered = None
        self.__taxonomy_df = taxonomy_df

    @property
    def rank_to_df(self):
        query_df = load_data.load_sample(self.metadata)
        if self.options.min_length:
            query_df = query_df[query_df['LENGTH'] >= self.options.min_length]
        rank_to_df = load_data.get_rank_to_df(query_df, self.taxonomy_df, self.label, self.is_gs)
        del query_df
        return rank_to_df

    @property
    def profile(self):
        if self.__profile:
            return self.__profile
        self.__profile = self._create_profile(all_bins=True)
        return self.__profile

    @property
    def profile_filtered(self):
        if self.__profile_filtered:
            return self.__profile_filtered
        self.__profile_filtered = self._create_profile(all_bins=False)
        return self.__profile_filtered

    @property
    def taxonomy_df(self):
        return self.__taxonomy_df

    @rank_to_df.setter
    def rank_to_df(self, rank_to_df):
        self.__rank_to_df = rank_to_df

    def _create_profile(self, all_bins):
        if self.recall_df.empty:
            return [], []

        class Prediction:
            def __init__(self):
                pass
        profile_bp = []
        profile_seq = []
        if all_bins:
            recall_df = self.recall_df
        else:
            recall_df = self.recall_df[~self.recall_df['filtered']]
        for index, row in recall_df.iterrows():
            prediction_bp = Prediction()
            prediction_bp.taxid = str(int(row['TAXID']))
            prediction_bp.rank = row['rank']
            prediction_bp.percentage = row['tp_length']
            taxpath = self.taxonomy_df.loc[row['TAXID']][load_ncbi_taxinfo.RANKS].values
            taxpath = '|'.join(map(lambda v: '' if pd.isnull(v) else str(v), taxpath)).rstrip('|')
            prediction_bp.taxpath = taxpath
            prediction_bp.taxpathsn = None
            profile_bp.append(prediction_bp)

            prediction_seq = Prediction()
            prediction_seq.taxid = prediction_bp.taxid
            prediction_seq.rank = prediction_bp.rank
            prediction_seq.percentage = row['tp_seq_counts']
            prediction_seq.taxpath = prediction_bp.taxpath
            prediction_seq.taxpathsn = None
            profile_seq.append(prediction_seq)
        return profile_bp, profile_seq

    def compute_unifrac(self, all_bins):
        if all_bins:
            pf_profile_bp = pf.Profile(profile=self.profile[0])
            pf_profile_seq = pf.Profile(profile=self.profile[1])
        else:
            pf_profile_bp = pf.Profile(profile=self.profile_filtered[0])
            pf_profile_seq = pf.Profile(profile=self.profile_filtered[1])
        gs_pf_profile_bp = pf.Profile(profile=self.gold_standard.profile[0])
        gs_pf_profile_seq = pf.Profile(profile=self.gold_standard.profile[1])
        return uf.compute_unifrac(gs_pf_profile_bp, pf_profile_bp), uf.compute_unifrac(gs_pf_profile_seq, pf_profile_seq)

    def get_metrics_df(self):
        allranks_metrics_df = pd.DataFrame()
        for rank in self.metrics:
            metrics_dict = self.metrics[rank].get_ordered_dict()
            rank_metrics_df = pd.DataFrame(metrics_dict)

            if self.metrics_filtered:
                rank_metrics_df = pd.DataFrame(metrics_dict) \
                    .drop(columns=[utils_labels.TOOL, utils_labels.BINNING_TYPE, utils_labels.SAMPLE, utils_labels.RANK]) \
                    .add_suffix(utils_labels.UNFILTERED)
                metrics_dict = self.metrics_filtered[rank].get_ordered_dict()
                rank_metrics_df_filtered = pd.DataFrame(metrics_dict)
                rank_metrics_df = pd.concat([rank_metrics_df_filtered, rank_metrics_df], axis=1)

            rank_metrics_df[utils_labels.TOOL] = self.label
            rank_metrics_df[utils_labels.BINNING_TYPE] = self.binning_type
            rank_metrics_df[utils_labels.RANK] = rank

            allranks_metrics_df = pd.concat([allranks_metrics_df, rank_metrics_df], ignore_index=True, sort=True)
        return allranks_metrics_df

    def compute_metrics_per_rank(self, rank_to_df, gs_rank_to_df, rank):
        if rank not in gs_rank_to_df:
            return
        logging.getLogger('amber').info('Evaluating {}, sample {}, taxonomic binning, {}'.format(self.label, self.sample_id, rank))

        self.metrics[rank] = Metrics()
        gs_df = gs_rank_to_df[rank]

        query_df = rank_to_df[rank][['SEQUENCEID', 'TAXID']]
        condition = query_df['SEQUENCEID'].isin(gs_df['SEQUENCEID'])
        if ~condition.all():
            logging.getLogger('amber').warning("{} sequences in {} not found in the gold standard.".format(query_df[~condition]['SEQUENCEID'].nunique(), self.label))
            query_df = query_df[condition]

        query_w_length_df = pd.merge(query_df, gs_df.rename(columns={'TAXID': 'true_taxid'}).reset_index(), on='SEQUENCEID', sort=False)

        confusion_df = query_w_length_df.groupby(['TAXID', 'true_taxid'], sort=False).agg({'LENGTH': 'sum', 'SEQUENCEID': 'count'}).rename(columns={'LENGTH': 'tax_length', 'SEQUENCEID': 'tax_seq_counts'})
        self.metrics[rank].rand_index_bp, self.metrics[rank].adjusted_rand_index_bp = Metrics.compute_rand_index(confusion_df, 'TAXID', 'true_taxid', 'tax_length')
        self.metrics[rank].rand_index_seq, self.metrics[rank].adjusted_rand_index_seq = Metrics.compute_rand_index(confusion_df, 'TAXID', 'true_taxid', 'tax_seq_counts')

        query_w_length_df = query_w_length_df[['SEQUENCEID', 'TAXID', 'LENGTH']]

        self.metrics[rank].percentage_of_assigned_bps = query_w_length_df['LENGTH'].sum() / gs_df['LENGTH'].sum()
        self.metrics[rank].percentage_of_assigned_seqs = query_w_length_df.shape[0] / gs_df.shape[0]

        true_positives_df = pd.merge(query_df, gs_df, on=['SEQUENCEID', 'TAXID'], sort=False)
        true_positives_df = true_positives_df.groupby('TAXID', sort=False).agg({'LENGTH': 'sum', 'SEQUENCEID': 'count'}).rename(columns={'LENGTH': 'tp_length', 'SEQUENCEID': 'tp_seq_counts'})

        tp_fp_fn_df = query_w_length_df.groupby('TAXID', sort=False).agg({'LENGTH': 'sum', 'SEQUENCEID': 'count'}).rename(columns={'LENGTH': 'total_length', 'SEQUENCEID': 'total_seq_counts'})
        tp_fp_fn_df = pd.merge(true_positives_df, tp_fp_fn_df, on=['TAXID'], how='outer', sort=False)
        tax_sizes_df = gs_df.groupby('TAXID', sort=False).agg({'LENGTH': 'sum', 'SEQUENCEID': 'count'}).rename(columns={'LENGTH': 'length_gs', 'SEQUENCEID': 'seq_counts_gs'})
        tp_fp_fn_df = tp_fp_fn_df.reset_index().join(tax_sizes_df, on='TAXID', how='outer', sort=False).set_index('TAXID').fillna(0).astype('int64')

        tp_fp_fn_df['precision_bp'] = tp_fp_fn_df['tp_length'] / tp_fp_fn_df['total_length']
        tp_fp_fn_df['precision_seq'] = tp_fp_fn_df['tp_seq_counts'] / tp_fp_fn_df['total_seq_counts']
        tp_fp_fn_df['recall_bp'] = tp_fp_fn_df['tp_length'] / tp_fp_fn_df['length_gs']
        tp_fp_fn_df['recall_seq'] = tp_fp_fn_df['tp_seq_counts'] / tp_fp_fn_df['seq_counts_gs']
        tp_fp_fn_df['rank'] = rank

        length_gs_sum = tp_fp_fn_df['length_gs'].sum()
        seq_counts_gs_sum = tp_fp_fn_df['seq_counts_gs'].sum()

        def set_values(metric_obj, df):
            tp_length_sum = df['tp_length'].sum()
            tp_seq_counts_sum = df['tp_seq_counts'].sum()

            metric_obj.precision_avg_bp = df['precision_bp'].mean()
            metric_obj.precision_avg_bp_sem = df['precision_bp'].sem()
            metric_obj.precision_avg_seq = df['precision_seq'].mean()
            metric_obj.precision_avg_seq_sem = df['precision_seq'].sem()
            metric_obj.precision_weighted_bp = tp_length_sum / df['total_length'].sum()
            metric_obj.precision_weighted_seq = tp_seq_counts_sum / df['total_seq_counts'].sum()

            metric_obj.recall_avg_bp = df['recall_bp'].mean()
            metric_obj.recall_avg_bp_sem = df['recall_bp'].sem()
            metric_obj.recall_avg_seq = df['recall_seq'].mean()
            metric_obj.recall_avg_seq_sem = df['recall_seq'].sem()
            metric_obj.recall_weighted_bp = tp_length_sum / length_gs_sum
            metric_obj.recall_weighted_seq = tp_seq_counts_sum / seq_counts_gs_sum

            metric_obj.accuracy_bp = tp_length_sum / length_gs_sum
            metric_obj.accuracy_seq = tp_seq_counts_sum / seq_counts_gs_sum

        set_values(self.metrics[rank], tp_fp_fn_df)

        if self.options.filter_tail_percentage:
            df_cpy = tp_fp_fn_df.copy()
            df_cpy['total_length_pct'] = df_cpy['total_length'] / df_cpy['total_length'].sum()
            df_cpy.sort_values(by='total_length', inplace=True)
            df_cpy['cumsum_length_pct'] = df_cpy['total_length_pct'].cumsum(axis=0)
            nan_rows = df_cpy['cumsum_length_pct'] <= self.options.filter_tail_percentage / 100
            df_cpy['precision_bp'].mask(nan_rows, inplace=True)
            df_cpy['precision_seq'].mask(nan_rows, inplace=True)
            df_cpy['recall_bp'].mask(nan_rows, other=.0, inplace=True)
            df_cpy['recall_seq'].mask(nan_rows, other=.0, inplace=True)
            df_cpy['tp_length'].mask(nan_rows, inplace=True)
            df_cpy['tp_seq_counts'].mask(nan_rows, inplace=True)
            df_cpy['total_length'].mask(nan_rows, inplace=True)
            df_cpy['total_seq_counts'].mask(nan_rows, inplace=True)
            df_cpy.drop(columns=['cumsum_length_pct', 'total_length_pct'], inplace=True)

            self.metrics_filtered[rank] = Metrics()
            set_values(self.metrics_filtered[rank], df_cpy)

            tp_fp_fn_df['filtered'] = nan_rows

            remaining = df_cpy[~nan_rows].index.values
            confusion_df = confusion_df.loc[remaining]
            self.metrics_filtered[rank].rand_index_bp, self.metrics_filtered[rank].adjusted_rand_index_bp = Metrics.compute_rand_index(confusion_df, 'TAXID', 'true_taxid', 'tax_length')
            self.metrics_filtered[rank].rand_index_seq, self.metrics_filtered[rank].adjusted_rand_index_seq = Metrics.compute_rand_index(confusion_df, 'TAXID', 'true_taxid', 'tax_seq_counts')

            query_w_length_df = query_w_length_df.set_index('TAXID').loc[remaining]
            self.metrics_filtered[rank].percentage_of_assigned_bps = query_w_length_df['LENGTH'].sum() / gs_df['LENGTH'].sum()
            self.metrics_filtered[rank].percentage_of_assigned_seqs = query_w_length_df.shape[0] / gs_df.shape[0]
        else:
            tp_fp_fn_df['filtered'] = False

        self.precision_df = pd.concat([self.precision_df, tp_fp_fn_df.reset_index().sort_values(
            by='recall_bp', axis=0, ascending=False)], ignore_index=True, sort=True)
        self.recall_df = self.precision_df

        self.recall_df['name'] = self.recall_df['TAXID'].apply(lambda x: self.taxonomy_df.loc[x]['name'])

    def compute_metrics(self, gs_rank_to_df):
        if self.label == utils_labels.GS and (self.options.only_genome_queries or self.options.skip_gs):
            return False
        rank_to_df = self.rank_to_df
        for rank in rank_to_df:
            self.compute_metrics_per_rank(rank_to_df, gs_rank_to_df, rank)
        del rank_to_df

        unifrac_bp, unifrac_seq = self.compute_unifrac(all_bins=True)
        for rank in self.metrics:
            self.metrics[rank].unifrac_bp = unifrac_bp
            self.metrics[rank].unifrac_seq = unifrac_seq
        if self.options.filter_tail_percentage:
            unifrac_bp, unifrac_seq = self.compute_unifrac(all_bins=False)
            for rank in self.metrics:
                self.metrics_filtered[rank].unifrac_bp = unifrac_bp
                self.metrics_filtered[rank].unifrac_seq = unifrac_seq

        return True


class Options:
    def __init__(self, filter_tail_percentage=0, genome_to_unique_common=None, filter_keyword=None, min_length=0,
                 rank_as_genome_binning=None, output_dir=None, min_completeness=None, max_contamination=None,
                 ncbi_dir=None, skip_gs=False):
        self.__filter_tail_percentage = float(filter_tail_percentage) if filter_tail_percentage else .0
        self.__genome_to_unique_common = genome_to_unique_common
        self.__filter_keyword = filter_keyword
        self.__min_length = int(min_length) if min_length else 0
        if rank_as_genome_binning and rank_as_genome_binning not in load_ncbi_taxinfo.RANKS:
            exit("Not a valid rank to assess taxonomic binning as genome binning (option --rank_as_genome_binning): " + rank_as_genome_binning)
        self.__rank_as_genome_binning = rank_as_genome_binning
        self.__only_genome_queries = True
        self.__only_taxonomic_queries = True
        self.__output_dir = output_dir
        if min_completeness:
            self.__min_completeness = [int(x.strip()) / 100.0 for x in min_completeness.split(',')]
        else:
            self.__min_completeness = [.5, .7, .9]
        if max_contamination:
            self.__max_contamination = [int(x.strip()) / 100.0 for x in max_contamination.split(',')]
        else:
            self.__max_contamination = [.1, .05]
        self.__skip_gs = skip_gs
        self.__ncbi_dir = ncbi_dir

    @property
    def filter_tail_percentage(self):
        return self.__filter_tail_percentage

    @property
    def genome_to_unique_common(self):
        return self.__genome_to_unique_common

    @property
    def filter_keyword(self):
        return self.__filter_keyword

    @property
    def min_length(self):
        return self.__min_length

    @property
    def rank_as_genome_binning(self):
        return self.__rank_as_genome_binning

    @property
    def only_genome_queries(self):
        return self.__only_genome_queries

    @property
    def only_taxonomic_queries(self):
        return self.__only_taxonomic_queries

    @property
    def output_dir(self):
        return self.__output_dir

    @property
    def min_completeness(self):
        return self.__min_completeness

    @property
    def max_contamination(self):
        return self.__max_contamination

    @property
    def ncbi_dir(self):
        return self.__ncbi_dir

    @property
    def skip_gs(self):
        return self.__skip_gs

    @filter_tail_percentage.setter
    def filter_tail_percentage(self, filter_tail_percentage):
        self.__filter_tail_percentage = filter_tail_percentage

    @genome_to_unique_common.setter
    def genome_to_unique_common(self, genome_to_unique_common):
        self.__genome_to_unique_common = genome_to_unique_common

    @filter_keyword.setter
    def filter_keyword(self, filter_keyword):
        self.__filter_keyword = filter_keyword

    @min_length.setter
    def min_length(self, min_length):
        self.__min_length = min_length

    @rank_as_genome_binning.setter
    def rank_as_genome_binning(self, rank_as_genome_binning):
        self.__rank_as_genome_binning = rank_as_genome_binning

    @only_genome_queries.setter
    def only_genome_queries(self, only_genome_queries):
        self.__only_genome_queries = only_genome_queries

    @only_taxonomic_queries.setter
    def only_taxonomic_queries(self, only_taxonomic_queries):
        self.__only_taxonomic_queries = only_taxonomic_queries

    @output_dir.setter
    def output_dir(self, output_dir):
        self.__output_dir = output_dir

    @min_completeness.setter
    def min_completeness(self, min_completeness):
        self.__min_completeness = min_completeness

    @max_contamination.setter
    def max_contamination(self, max_contamination):
        self.__max_contamination = max_contamination

    @ncbi_dir.setter
    def ncbi_dir(self, ncbi_dir):
        self.__ncbi_dir = ncbi_dir

    @skip_gs.setter
    def skip_gs(self, skip_gs):
        self.__skip_gs = skip_gs
