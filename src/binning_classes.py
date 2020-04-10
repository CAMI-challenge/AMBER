# Copyright 2020 Department of Computational Biology for Infection Research - Helmholtz Centre for Infection Research
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
from abc import ABC, abstractmethod
from collections import defaultdict
from collections import OrderedDict
from src.utils import labels as utils_labels


class Metrics():
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
        self.__precision_avg_seq  = .0
        self.__precision_avg_seq_sem = .0
        self.__precision_weighted_bp = .0
        self.__precision_weighted_seq = .0
        self.__recall_avg_bp = .0
        self.__recall_avg_bp_var = .0
        self.__recall_avg_bp_sem = .0
        self.__recall_avg_seq = .0
        self.__recall_avg_seq_sem = .0
        self.__recall_weighted_bp = .0
        self.__recall_weighted_seq = .0
        self.__f1_score_bp = .0
        self.__f1_score_seq = .0

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
    def recall_avg_bp_var(self):
        return self.__recall_avg_bp_var

    @property
    def recall_avg_bp_sem(self):
        return self.__recall_avg_bp_sem

    @property
    def recall_avg_seq(self):
        return self.__recall_avg_seq

    @property
    def recall_avg_seq_sem(self):
        return self.__recall_avg_seq_sem

    @property
    def recall_weighted_bp(self):
        return self.__recall_weighted_bp

    @property
    def recall_weighted_seq(self):
        return self.__recall_weighted_seq

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

    @recall_avg_bp_var.setter
    def recall_avg_bp_var(self, recall_avg_bp_var):
        self.__recall_avg_bp_var = recall_avg_bp_var

    @recall_avg_bp_sem.setter
    def recall_avg_bp_sem(self, recall_avg_bp_sem):
        self.__recall_avg_bp_sem = recall_avg_bp_sem

    @recall_avg_seq.setter
    def recall_avg_seq(self, recall_avg_seq):
        self.__recall_avg_seq = recall_avg_seq

    @recall_avg_seq_sem.setter
    def recall_avg_seq_sem(self, recall_avg_seq_sem):
        self.__recall_avg_seq_sem = recall_avg_seq_sem

    @recall_weighted_bp.setter
    def recall_weighted_bp(self, recall_weighted_bp):
        self.__recall_weighted_bp = recall_weighted_bp

    @recall_weighted_seq.setter
    def recall_weighted_seq(self, recall_weighted_seq):
        self.__recall_weighted_seq = recall_weighted_seq

    @staticmethod
    def compute_rand_index(confusion_df, col_name, gs_col_name, field):
        def choose2(n):
            return (n * (n - 1)) / 2.0

        bin_mapping_comb = confusion_df[field].apply(choose2).sum()
        bin_comb = confusion_df.groupby(col_name).agg({field: 'sum'})[field].apply(choose2).sum()
        mapping_comb = confusion_df.groupby(gs_col_name).agg({field: 'sum'})[field].apply(choose2).sum()
        num_bp_comb = choose2(confusion_df[field].sum())

        rand_index = ((num_bp_comb - bin_comb - mapping_comb + 2 * bin_mapping_comb) / num_bp_comb) if num_bp_comb != 0 else .0

        temp = (bin_comb * mapping_comb / num_bp_comb) if num_bp_comb != 0 else .0
        ret = bin_mapping_comb - temp
        denominator = (((bin_comb + mapping_comb) / 2.0) - temp)
        adjusted_rand_index = (ret / denominator) if denominator != 0 else .0

        return rand_index, adjusted_rand_index

    def get_ordered_dict(self):
        return OrderedDict([(utils_labels.TOOL, None),
                            (utils_labels.BINNING_TYPE, None),
                            (utils_labels.SAMPLE, None),
                            (utils_labels.RANK, None),
                            (utils_labels.AVG_PRECISION_BP, [self.__precision_avg_bp]),
                            (utils_labels.AVG_PRECISION_BP_SEM, [self.__precision_avg_bp_sem]),

                            ('avg_precision_bp_var', [self.__precision_avg_bp_var]),
                            (utils_labels.AVG_RECALL_BP, [self.__recall_avg_bp]),
                            (utils_labels.AVG_RECALL_BP_SEM, [self.__recall_avg_bp_sem]),
                            ('avg_recall_bp_var', [self.__recall_avg_bp_var]),
                            (utils_labels.F1_SCORE_BP, [2 * self.__precision_avg_bp * self.__recall_avg_bp / (self.__precision_avg_bp + self.__recall_avg_bp)]),

                            (utils_labels.AVG_PRECISION_SEQ, [self.__precision_avg_seq]),
                            (utils_labels.AVG_PRECISION_SEQ_SEM, [self.__precision_avg_seq_sem]),
                            (utils_labels.AVG_RECALL_SEQ, [self.__recall_avg_seq]),
                            (utils_labels.AVG_RECALL_SEQ_SEM, [self.__recall_avg_seq_sem]),
                            (utils_labels.F1_SCORE_SEQ, [2 * self.__precision_avg_seq * self.__recall_avg_seq / (self.__precision_avg_seq + self.__recall_avg_seq)]),

                            (utils_labels.PRECISION_PER_BP, [self.__precision_weighted_bp]),
                            (utils_labels.PRECISION_PER_SEQ, [self.__precision_weighted_seq]),
                            (utils_labels.RECALL_PER_BP, [self.__recall_weighted_bp]),
                            (utils_labels.RECALL_PER_SEQ, [self.__recall_weighted_seq]),
                            (utils_labels.F1_SCORE_PER_BP, [2 * self.__precision_weighted_bp * self.__recall_weighted_bp / (self.__precision_weighted_bp + self.__recall_weighted_bp)]),
                            (utils_labels.F1_SCORE_PER_SEQ, [2 * self.__precision_weighted_seq * self.__recall_weighted_seq / (self.__precision_weighted_seq + self.__recall_weighted_seq)]),

                            (utils_labels.ACCURACY_PER_BP, [self.__accuracy_bp]),
                            (utils_labels.ACCURACY_PER_SEQ, [self.__accuracy_seq]),

                            (utils_labels.PERCENTAGE_ASSIGNED_BPS, [self.__percentage_of_assigned_bps]),
                            (utils_labels.PERCENTAGE_ASSIGNED_SEQS, [self.__percentage_of_assigned_seqs]),
                            (utils_labels.RI_BY_BP, [self.__rand_index_bp]),
                            (utils_labels.RI_BY_SEQ, [self.__rand_index_seq]),
                            (utils_labels.ARI_BY_BP, [self.__adjusted_rand_index_bp]),
                            (utils_labels.ARI_BY_SEQ, [self.__adjusted_rand_index_seq]),

                            (utils_labels.UNIFRAC_BP, [0]),
                            (utils_labels.UNIFRAC_SEQ, [0]),

                            (utils_labels.MISCLASSIFICATION_PER_BP, [1 - self.__precision_weighted_bp]),
                            (utils_labels.MISCLASSIFICATION_PER_SEQ, [1 - self.__precision_weighted_seq])])


class Query(ABC):
    def __init__(self, label):
        self.__label = label
        self.__gold_standard_df = None
        self.__precision_df = pd.DataFrame()
        self.__recall_df = None
        self.__confusion_df = None
        self.__metrics = None

    @property
    def label(self):
        return self.__label

    @property
    def gold_standard_df(self):
        return self.__gold_standard_df

    @property
    def precision_df(self):
        return self.__precision_df

    @property
    def recall_df(self):
        return self.__recall_df

    @property
    def confusion_df(self):
        return self.__confusion_df

    @property
    def metrics(self):
        return self.__metrics

    @label.setter
    def label(self, label):
        self.__label = label

    @gold_standard_df.setter
    def gold_standard_df(self, gold_standard_df):
        self.__gold_standard_df = gold_standard_df

    @precision_df.setter
    def precision_df(self, precision_df):
        self.__precision_df = precision_df

    @recall_df.setter
    def recall_df(self, recall_df):
        self.__recall_df = recall_df

    @confusion_df.setter
    def confusion_df(self, confusion_df):
        self.__confusion_df = confusion_df

    @metrics.setter
    def metrics(self, metrics):
        self.__metrics = metrics

    @abstractmethod
    def compute_metrics(self):
        pass


class GenomeQuery(Query):
    binning_type = 'genome'

    def __init__(self, df, label):
        super().__init__(label)
        self.__df = df
        self.metrics = Metrics()

    @property
    def df(self):
        return self.__df

    @df.setter
    def df(self, df):
        self.__df = df

    def get_metrics_df(self):
        metrics_dict = self.metrics.get_ordered_dict()
        metrics_dict[utils_labels.TOOL] = self.label
        metrics_dict[utils_labels.BINNING_TYPE] = self.binning_type
        metrics_dict[utils_labels.RANK] = 'NA'
        return pd.DataFrame(metrics_dict)

    def compute_metrics(self):
        print(self.label)
        gs_df = self.gold_standard_df
        # gs_df.rename(columns={'BINID': 'bin_id', 'LENGTH': 'seq_length'}, inplace=True)
        gs_df = gs_df[['SEQUENCEID', 'BINID', 'LENGTH']].rename(columns={'LENGTH': 'seq_length'})

        query_df = self.df
        # query_df.rename(columns={'BINID': 'bin_id'}, inplace=True)
        query_df = query_df[['SEQUENCEID', 'BINID']]

        query_w_length = pd.merge(query_df, gs_df[['SEQUENCEID', 'BINID', 'seq_length']].rename(columns={'BINID': 'genome_id'}).drop_duplicates('SEQUENCEID'), on='SEQUENCEID', sort=False)

        query_w_length_no_dups = query_w_length.drop_duplicates('SEQUENCEID')
        gs_df_no_dups = gs_df.drop_duplicates('SEQUENCEID')
        self.metrics.percentage_of_assigned_bps = query_w_length_no_dups['seq_length'].sum() / gs_df_no_dups['seq_length'].sum()
        self.metrics.percentage_of_assigned_seqs = query_w_length_no_dups.shape[0] / gs_df_no_dups['SEQUENCEID'].shape[0]

        # confusion table possibly with the same sequences in multiple bins
        query_w_length_mult_seqs = query_df.reset_index().merge(gs_df, on='SEQUENCEID', sort=False) #.set_index('index') #.set_index(['index', 'SEQUENCEID'])

        if query_w_length.shape[0] < query_w_length_mult_seqs.shape[0]:
            query_w_length_mult_seqs.drop_duplicates(['index', 'genome_id'], inplace=True)
            confusion_df = query_w_length_mult_seqs.groupby(['BINID', 'genome_id'], sort=False).agg({'seq_length': 'sum', 'SEQUENCEID': 'count'}).rename(columns={'seq_length': 'genome_length', 'SEQUENCEID': 'genome_seq_counts'})

            most_abundant_genome_df = confusion_df.loc[confusion_df.groupby('BINID', sort=False)['genome_length'].idxmax()]
            most_abundant_genome_df = most_abundant_genome_df.reset_index()[['BINID', 'genome_id']]

            matching_genomes_df = pd.merge(query_w_length_mult_seqs, most_abundant_genome_df, on=['BINID', 'genome_id']).set_index('index')
            query_w_length_mult_seqs.set_index('index', inplace=True)
            difference_df = query_w_length_mult_seqs.drop(matching_genomes_df.index).groupby(['index'], sort=False).first()
            query_w_length = pd.concat([matching_genomes_df, difference_df])

            # query_w_length_mult_seqs.reset_index(inplace=True)
            # query_w_length_mult_seqs = pd.merge(query_w_length_mult_seqs, most_abundant_genome_df, on=['BINID'])
            # grouped = query_w_length_mult_seqs.groupby(['index'], sort=False, as_index=False)
            # query_w_length = grouped.apply(lambda x: x[x['genome_id_x'] == x['genome_id_y'] if any(x['genome_id_x'] == x['genome_id_y']) else len(x) * [True]])
            # query_w_length = query_w_length.groupby(['index'], sort=False).first().drop(columns='genome_id_y').rename(columns={'genome_id_x': 'genome_id'})

        self.df = query_w_length

        confusion_df = query_w_length.groupby(['BINID', 'genome_id'], sort=False).agg({'seq_length': 'sum', 'SEQUENCEID': 'count'}).rename(columns={'seq_length': 'genome_length', 'SEQUENCEID': 'genome_seq_counts'})
        self.confusion_df = confusion_df

        self.metrics.rand_index_bp, self.metrics.adjusted_rand_index_bp = Metrics.compute_rand_index(confusion_df, 'BINID', 'genome_id', 'genome_length')
        self.metrics.rand_index_seq, self.metrics.adjusted_rand_index_seq = Metrics.compute_rand_index(confusion_df, 'BINID', 'genome_id', 'genome_seq_counts')

        most_abundant_genome_df = confusion_df.loc[confusion_df.groupby('BINID', sort=False)['genome_length'].idxmax()].reset_index().set_index('BINID')

        precision_df = query_w_length.groupby('BINID', sort=False).agg({'seq_length': 'sum', 'SEQUENCEID': 'count'}).rename(columns={'seq_length': 'total_length', 'SEQUENCEID': 'total_seq_counts'})
        precision_df = pd.merge(precision_df, most_abundant_genome_df, on='BINID')
        precision_df['precision_bp'] = precision_df['genome_length'] / precision_df['total_length']
        precision_df['precision_seq'] = precision_df['genome_seq_counts'] / precision_df['total_seq_counts']

        self.metrics.precision_avg_bp = precision_df['precision_bp'].mean()
        self.metrics.precision_avg_bp_sem = precision_df['precision_bp'].sem()
        self.metrics.precision_avg_bp_var = precision_df['precision_bp'].var()
        self.metrics.precision_avg_seq = precision_df['precision_seq'].mean()
        self.metrics.precision_weighted_bp = precision_df['genome_length'].sum() / precision_df['total_length'].sum()
        self.metrics.precision_weighted_seq = precision_df['genome_seq_counts'].sum() / precision_df['total_seq_counts'].sum()

        genome_sizes_df = gs_df.rename(columns={'BINID': 'genome_id'}).groupby('genome_id', sort=False).agg({'seq_length': 'sum', 'SEQUENCEID': 'count'}).rename(columns={'seq_length': 'length_gs', 'SEQUENCEID': 'seq_counts_gs'})
        precision_df = precision_df.reset_index().join(genome_sizes_df, on='genome_id', how='left', sort=False).set_index('BINID')
        precision_df['recall_bp'] = precision_df['genome_length'] / precision_df['length_gs']
        precision_df['recall_seq'] = precision_df['genome_seq_counts'] / precision_df['seq_counts_gs']
        precision_df['rank'] = 'NA'

        recall_df = confusion_df.loc[confusion_df.groupby('genome_id', sort=False)['genome_length'].idxmax()]
        recall_df = recall_df.reset_index().join(genome_sizes_df, on='genome_id', how='right', sort=False).set_index('BINID')
        recall_df.fillna({'genome_length': 0, 'genome_seq_counts': 0}, inplace=True)
        recall_df['recall_bp'] = recall_df['genome_length'] / recall_df['length_gs']
        recall_df['recall_seq'] = recall_df['genome_seq_counts'] / recall_df['seq_counts_gs']

        self.metrics.recall_avg_bp = recall_df['recall_bp'].mean()
        self.metrics.recall_avg_bp_var = recall_df['recall_bp'].var()
        self.metrics.recall_avg_bp_sem = recall_df['recall_bp'].sem()
        self.metrics.recall_avg_seq = recall_df['recall_seq'].mean()
        self.metrics.recall_avg_seq_sem = recall_df['recall_seq'].sem()
        self.metrics.recall_weighted_bp = recall_df['genome_length'].sum() / recall_df['length_gs'].sum()
        self.metrics.recall_weighted_seq = recall_df['genome_seq_counts'].sum() / recall_df['seq_counts_gs'].sum()

        self.metrics.accuracy_bp = precision_df['genome_length'].sum() / recall_df['length_gs'].sum()
        self.metrics.accuracy_seq = precision_df['genome_seq_counts'].sum() / recall_df['seq_counts_gs'].sum()

        self.precision_df = precision_df
        self.recall_df = recall_df


class TaxonomicQuery(Query):
    tax_id_to_parent = None
    tax_id_to_rank = None
    tax_id_to_name = None
    tax_id_to_tax_id = None
    binning_type = 'taxonomic'

    def __init__(self, rank_to_df, label):
        super().__init__(label)
        self.__rank_to_df = rank_to_df
        self.metrics = defaultdict()


    @property
    def rank_to_df(self):
        return self.__rank_to_df

    @rank_to_df.setter
    def rank_to_df(self, rank_to_df):
        self.__rank_to_df = rank_to_df

    def get_metrics_df(self):
        allranks_metrics_df = pd.DataFrame()
        for rank in self.metrics:
            metrics_dict = self.metrics[rank].get_ordered_dict()
            metrics_dict[utils_labels.TOOL] = self.label
            metrics_dict[utils_labels.BINNING_TYPE] = self.binning_type
            metrics_dict[utils_labels.RANK] = rank
            rank_metrics_df  = pd.DataFrame(metrics_dict)
            allranks_metrics_df = pd.concat([allranks_metrics_df, rank_metrics_df], ignore_index=True, sort=True)
        return allranks_metrics_df

    def compute_metrics_for_rank(self, rank):
        self.metrics[rank] = Metrics()
        gs_df = self.gold_standard_df[rank].reset_index()

        query_df = self.rank_to_df[rank].reset_index()[['SEQUENCEID', 'TAXID']]

        query_w_length_df = pd.merge(query_df, gs_df.rename(columns={'TAXID': 'true_taxid'}).reset_index(),  on='SEQUENCEID', sort=False)

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
        tp_fp_fn_df = tp_fp_fn_df.reset_index().join(tax_sizes_df, on='TAXID', how='outer', sort=False).set_index('TAXID')
        tp_fp_fn_df.fillna(0, inplace=True)

        tp_fp_fn_df['precision_bp'] = tp_fp_fn_df['tp_length'] / tp_fp_fn_df['total_length']
        tp_fp_fn_df['precision_seq'] = tp_fp_fn_df['tp_seq_counts'] / tp_fp_fn_df['total_seq_counts']
        tp_fp_fn_df['recall_bp'] = tp_fp_fn_df['tp_length'] / tp_fp_fn_df['length_gs']
        tp_fp_fn_df['recall_seq'] = tp_fp_fn_df['tp_seq_counts'] / tp_fp_fn_df['seq_counts_gs']
        tp_fp_fn_df['rank'] = rank

        tp_length_sum = tp_fp_fn_df['tp_length'].sum()
        tp_seq_counts_sum = tp_fp_fn_df['tp_seq_counts'].sum()
        length_gs_sum = tp_fp_fn_df['length_gs'].sum()
        seq_counts_gs_sum = tp_fp_fn_df['seq_counts_gs'].sum()

        self.metrics[rank].precision_avg_bp = tp_fp_fn_df['precision_bp'].mean()
        self.metrics[rank].precision_avg_seq = tp_fp_fn_df['precision_seq'].mean()
        self.metrics[rank].precision_weighted_bp = tp_length_sum / tp_fp_fn_df['total_length'].sum()
        self.metrics[rank].precision_weighted_seq = tp_seq_counts_sum / tp_fp_fn_df['total_seq_counts'].sum()

        self.metrics[rank].recall_avg_bp = tp_fp_fn_df['recall_bp'].mean()
        self.metrics[rank].recall_avg_seq = tp_fp_fn_df['recall_seq'].mean()
        self.metrics[rank].recall_weighted_bp = tp_length_sum / length_gs_sum
        self.metrics[rank].recall_weighted_seq = tp_seq_counts_sum / seq_counts_gs_sum

        self.metrics[rank].accuracy_bp = tp_length_sum / length_gs_sum
        self.metrics[rank].accuracy_seq = tp_seq_counts_sum / seq_counts_gs_sum

        self.precision_df = pd.concat([self.precision_df, tp_fp_fn_df.reset_index()], ignore_index=True, sort=True)
        self.recall_df = self.precision_df

    def compute_metrics(self):
        print(self.label)
        for rank in self.rank_to_df:
            self.compute_metrics_for_rank(rank)
