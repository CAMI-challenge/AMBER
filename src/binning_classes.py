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
import logging
import itertools
from abc import ABC, abstractmethod
from collections import defaultdict
from collections import OrderedDict
from src.utils import labels as utils_labels
from src.utils import load_ncbi_taxinfo
from src.utils import ProfilingTools as pf
from src import unifrac_distance as uf


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

    @unifrac_bp.setter
    def unifrac_bp(self, unifrac_bp):
        self.__unifrac_bp = unifrac_bp

    @unifrac_seq.setter
    def unifrac_seq(self, unifrac_seq):
        self.__unifrac_seq = unifrac_seq

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

                            (utils_labels.UNIFRAC_BP, self.__unifrac_bp),
                            (utils_labels.UNIFRAC_SEQ, self.__unifrac_seq),

                            (utils_labels.MISCLASSIFICATION_PER_BP, [1 - self.__precision_weighted_bp]),
                            (utils_labels.MISCLASSIFICATION_PER_SEQ, [1 - self.__precision_weighted_seq])])


class Query(ABC):
    def __init__(self, label, options):
        self.__label = label
        self.__gold_standard = None
        self.__gold_standard_df = None
        self.__precision_df = pd.DataFrame()
        self.__recall_df = None
        self.__confusion_df = None
        self.__metrics = None
        self.__options = options

    @property
    def label(self):
        return self.__label

    @property
    def gold_standard(self):
        return self.__gold_standard

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

    @property
    def options(self):
        return self.__options

    @label.setter
    def label(self, label):
        self.__label = label

    @gold_standard.setter
    def gold_standard(self, gold_standard):
        self.__gold_standard = gold_standard

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

    @options.setter
    def options(self, options):
        self.__options = options

    @abstractmethod
    def compute_metrics(self):
        pass

    def compute_unifrac(self):
        return None, None


class GenomeQuery(Query):
    binning_type = 'genome'

    def __init__(self, df, label, options):
        super().__init__(label, options)
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
        if self.label == utils_labels.GS and self.options.only_taxonomic_queries:
            return False
        logging.getLogger('amber').info('Evaluating {} (genome binning)'.format(self.label))
        gs_df = self.gold_standard_df
        gs_df = gs_df[['SEQUENCEID', 'BINID', 'LENGTH']].rename(columns={'LENGTH': 'seq_length'})

        query_df = self.df
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

        if self.options.genome_to_unique_common:
            recall_df = recall_df[~recall_df['genome_id'].isin(self.options.genome_to_unique_common)]

        self.metrics.recall_avg_bp = recall_df['recall_bp'].mean()
        self.metrics.recall_avg_bp_var = recall_df['recall_bp'].var()
        self.metrics.recall_avg_bp_sem = recall_df['recall_bp'].sem()
        self.metrics.recall_avg_seq = recall_df['recall_seq'].mean()
        self.metrics.recall_avg_seq_sem = recall_df['recall_seq'].sem()
        self.metrics.recall_weighted_bp = recall_df['genome_length'].sum() / recall_df['length_gs'].sum()
        self.metrics.recall_weighted_seq = recall_df['genome_seq_counts'].sum() / recall_df['seq_counts_gs'].sum()

        self.metrics.accuracy_bp = precision_df['genome_length'].sum() / recall_df['length_gs'].sum()
        self.metrics.accuracy_seq = precision_df['genome_seq_counts'].sum() / recall_df['seq_counts_gs'].sum()

        self.precision_df = precision_df.sort_values(by=['recall_bp'], axis=0, ascending=False)
        self.recall_df = recall_df
        return True

    @staticmethod
    def calc_num_recovered_genomes(pd_bins, min_completeness, max_contamination):
        counts_list = []
        for (sample_id, tool), pd_group in pd_bins.groupby(['sample_id', utils_labels.TOOL]):
            for x in itertools.product(min_completeness, max_contamination):
                count = pd_group[(pd_group['recall_bp'] > x[0]) & (pd_group['precision_bp'] > (1 - x[1]))].shape[0]
                counts_list.append((sample_id, tool, '> ' + str(x[0]) + '% completeness', '< ' + str(x[1]) + '%', count))

        pd_counts = pd.DataFrame(counts_list, columns=[utils_labels.SAMPLE, utils_labels.TOOL, 'Completeness', 'Contamination', 'count'])
        pd_counts = pd.pivot_table(pd_counts, values='count', index=[utils_labels.SAMPLE, utils_labels.TOOL, 'Contamination'], columns=['Completeness']).reset_index()
        return pd_counts


class TaxonomicQuery(Query):
    tax_id_to_parent = None
    tax_id_to_rank = None
    tax_id_to_name = None
    tax_id_to_tax_id = None
    binning_type = 'taxonomic'

    def __init__(self, rank_to_df, label, options):
        super().__init__(label, options)
        self.__rank_to_df = rank_to_df
        self.metrics = defaultdict()
        self.__profile = None

    @property
    def rank_to_df(self):
        return self.__rank_to_df

    @property
    def profile(self):
        if self.__profile:
            return self.__profile
        self.__profile = self._create_profile()
        return self.__profile

    @rank_to_df.setter
    def rank_to_df(self, rank_to_df):
        self.__rank_to_df = rank_to_df

    def _create_profile(self):
        class Prediction:
            def __init__(self):
                pass
        profile_bp = []
        profile_seq = []
        for index, row in self.recall_df.iterrows():
            prediction_bp = Prediction()
            prediction_bp.taxid = str(row['TAXID'])
            prediction_bp.rank = row['rank']
            prediction_bp.percentage = row['tp_length']
            taxpath = load_ncbi_taxinfo.get_id_path(row['TAXID'], TaxonomicQuery.tax_id_to_parent, TaxonomicQuery.tax_id_to_rank, None)
            prediction_bp.taxpath = '|'.join(map(str, taxpath))
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

    def compute_unifrac(self):
        pf_profile_bp = pf.Profile(profile=self.profile[0])
        gs_pf_profile_bp = pf.Profile(profile=self.gold_standard.profile[0])
        pf_profile_seq = pf.Profile(profile=self.profile[1])
        gs_pf_profile_seq = pf.Profile(profile=self.gold_standard.profile[1])
        return uf.compute_unifrac(gs_pf_profile_bp, pf_profile_bp), uf.compute_unifrac(gs_pf_profile_seq, pf_profile_seq)

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

        if self.options.filter_tail_percentage:
            tp_fp_fn_df['total_length_pct'] = tp_fp_fn_df['total_length'] / tp_fp_fn_df['total_length'].sum()
            tp_fp_fn_df.sort_values(by='total_length', inplace=True)
            tp_fp_fn_df['cumsum_length_pct'] = tp_fp_fn_df['total_length_pct'].cumsum(axis=0)
            tp_fp_fn_df['precision_bp'].mask(tp_fp_fn_df['cumsum_length_pct'] <= self.options.filter_tail_percentage / 100, inplace=True)
            tp_fp_fn_df['precision_seq'].mask(tp_fp_fn_df['precision_bp'].isna(), inplace=True)
            tp_fp_fn_df.drop(columns=['cumsum_length_pct', 'total_length_pct'], inplace=True)

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

        self.precision_df = pd.concat([self.precision_df, tp_fp_fn_df.reset_index().sort_values(
            by='recall_bp', axis=0, ascending=False)], ignore_index=True, sort=True)
        self.recall_df = self.precision_df

        if self.tax_id_to_name:
            self.recall_df['name'] = self.recall_df['TAXID'].map(self.tax_id_to_name)

    def compute_metrics(self):
        if self.label == utils_labels.GS and self.options.only_genome_queries:
            return False
        logging.getLogger('amber').info('Evaluating {} (taxonomic binning)'.format(self.label))
        for rank in self.rank_to_df:
            self.compute_metrics_for_rank(rank)

        unifrac_bp, unifrac_seq = self.compute_unifrac()
        for rank in self.metrics:
            self.metrics[rank].unifrac_bp = unifrac_bp
            self.metrics[rank].unifrac_seq = unifrac_seq
        return True


class Options:
    def __init__(self, filter_tail_percentage, genome_to_unique_common, filter_keyword, min_length,
                 rank_as_genome_binning, output_dir, min_completeness=None, max_contamination=None):
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
