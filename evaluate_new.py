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
from collections import OrderedDict
from src import new_classes
from src.utils import labels as utils_labels


def get_rand_index(confusion_df, col_name, gs_col_name, field):
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


def compute_metrics(query):
    print(query.label)
    gs_df = query.gold_standard_df
    gs_df.rename(columns={'BINID': 'bin_id', 'LENGTH': 'seq_length'}, inplace=True)
    gs_df = gs_df[['SEQUENCEID', 'bin_id', 'seq_length']]

    query_df = query.df
    query_df.rename(columns={'BINID': 'bin_id'}, inplace=True)
    query_df = query_df[['SEQUENCEID', 'bin_id']]

    query_w_length = pd.merge(query_df, gs_df[['SEQUENCEID', 'bin_id', 'seq_length']].rename(columns={'bin_id': 'genome_id'}).drop_duplicates('SEQUENCEID'), on='SEQUENCEID', sort=False)

    query_w_length_no_dups = query_w_length.drop_duplicates('SEQUENCEID')
    gs_df_no_dups = gs_df.drop_duplicates('SEQUENCEID')
    percentage_of_assigned_bps = query_w_length_no_dups['seq_length'].sum() / gs_df_no_dups['seq_length'].sum()
    percentage_of_assigned_seqs = query_w_length_no_dups.shape[0] / gs_df_no_dups['SEQUENCEID'].shape[0]

    # confusion table possibly with the same sequences in multiple bins
    query_w_length_mult_seqs = query_df.reset_index().merge(gs_df, on='SEQUENCEID', sort=False) #.set_index('index') #.set_index(['index', 'SEQUENCEID'])

    if query_w_length.shape[0] < query_w_length_mult_seqs.shape[0]:
        query_w_length_mult_seqs.drop_duplicates(['index', 'genome_id'], inplace=True)
        confusion_df = query_w_length_mult_seqs.groupby(['bin_id', 'genome_id'], sort=False).agg({'seq_length': 'sum', 'SEQUENCEID': 'count'}).rename(columns={'seq_length': 'genome_length', 'SEQUENCEID': 'genome_seq_counts'})

        most_abundant_genome_df = confusion_df.loc[confusion_df.groupby('bin_id', sort=False)['genome_length'].idxmax()]
        most_abundant_genome_df = most_abundant_genome_df.reset_index()[['bin_id', 'genome_id']]

        matching_genomes_df = pd.merge(query_w_length_mult_seqs, most_abundant_genome_df, on=['bin_id', 'genome_id']).set_index('index')
        query_w_length_mult_seqs.set_index('index', inplace=True)
        difference_df = query_w_length_mult_seqs.drop(matching_genomes_df.index).groupby(['index'], sort=False).first()
        query_w_length = pd.concat([matching_genomes_df, difference_df])

        # query_w_length_mult_seqs.reset_index(inplace=True)
        # query_w_length_mult_seqs = pd.merge(query_w_length_mult_seqs, most_abundant_genome_df, on=['bin_id'])
        # grouped = query_w_length_mult_seqs.groupby(['index'], sort=False, as_index=False)
        # query_w_length = grouped.apply(lambda x: x[x['genome_id_x'] == x['genome_id_y'] if any(x['genome_id_x'] == x['genome_id_y']) else len(x) * [True]])
        # query_w_length = query_w_length.groupby(['index'], sort=False).first().drop(columns='genome_id_y').rename(columns={'genome_id_x': 'genome_id'})

    query.df = query_w_length

    confusion_df = query_w_length.groupby(['bin_id', 'genome_id'], sort=False).agg({'seq_length': 'sum', 'SEQUENCEID': 'count'}).rename(columns={'seq_length': 'genome_length', 'SEQUENCEID': 'genome_seq_counts'})
    query.confusion_df = confusion_df

    rand_index_bp, adjusted_rand_index_bp = get_rand_index(confusion_df, 'bin_id', 'genome_id', 'genome_length')
    rand_index_seq, adjusted_rand_index_seq = get_rand_index(confusion_df, 'bin_id', 'genome_id', 'genome_seq_counts')

    most_abundant_genome_df = confusion_df.loc[confusion_df.groupby('bin_id', sort=False)['genome_length'].idxmax()].reset_index().set_index('bin_id')

    precision_df = query_w_length.groupby('bin_id', sort=False).agg({'seq_length': 'sum', 'SEQUENCEID': 'count'}).rename(columns={'seq_length': 'total_length', 'SEQUENCEID': 'total_seq_counts'})
    # print(precision_df)
    precision_df = pd.merge(precision_df, most_abundant_genome_df, on='bin_id')
    # print(precision_df)
    precision_df['precision_bp'] = precision_df['genome_length'] / precision_df['total_length']
    precision_df['precision_seq'] = precision_df['genome_seq_counts'] / precision_df['total_seq_counts']
    # exit()
    precision_avg_bp = precision_df['precision_bp'].mean()
    precision_avg_seq = precision_df['precision_seq'].mean()
    precision_weighted_bp = precision_df['genome_length'].sum() / precision_df['total_length'].sum()
    precision_weighted_seq = precision_df['genome_seq_counts'].sum() / precision_df['total_seq_counts'].sum()


    genome_sizes_df = gs_df.rename(columns={'bin_id': 'genome_id'}).groupby('genome_id', sort=False).agg({'seq_length': 'sum', 'SEQUENCEID': 'count'}).rename(columns={'seq_length': 'length_gs', 'SEQUENCEID': 'seq_counts_gs'})
    precision_df = precision_df.reset_index().join(genome_sizes_df, on='genome_id', how='left', sort=False).set_index('bin_id')
    precision_df['recall_bp'] = precision_df['genome_length'] / precision_df['length_gs']
    precision_df['recall_seq'] = precision_df['genome_seq_counts'] / precision_df['seq_counts_gs']


    recall_df = confusion_df.loc[confusion_df.groupby('genome_id', sort=False)['genome_length'].idxmax()]
    recall_df = recall_df.reset_index().join(genome_sizes_df, on='genome_id', how='right', sort=False).set_index('bin_id')
    recall_df.fillna({'genome_length': 0, 'genome_seq_counts': 0}, inplace=True)
    recall_df['recall_bp'] = recall_df['genome_length'] / recall_df['length_gs']
    recall_df['recall_seq'] = recall_df['genome_seq_counts'] / recall_df['seq_counts_gs']

    recall_avg_bp = recall_df['recall_bp'].mean()
    recall_avg_seq = recall_df['recall_seq'].mean()
    recall_weighted_bp = recall_df['genome_length'].sum() / recall_df['length_gs'].sum()
    recall_weighted_seq = recall_df['genome_seq_counts'].sum() / recall_df['seq_counts_gs'].sum()

    accuracy_bp = precision_df['genome_length'].sum() / recall_df['length_gs'].sum()
    accuracy_seq = precision_df['genome_seq_counts'].sum() / recall_df['seq_counts_gs'].sum()

    query.precision_df = precision_df
    query.recall_df = recall_df

    f1_score_bp = 2 * precision_avg_bp * recall_avg_bp / (precision_avg_bp + recall_avg_bp)
    # f1_score_weighted_bp = 2 * precision_weighted_bp * recall_weighted_bp / (precision_weighted_bp + recall_weighted_bp)
    print('f1_score_bp: {}'.format(f1_score_bp))

    return percentage_of_assigned_bps, percentage_of_assigned_seqs, accuracy_bp, accuracy_seq, \
           rand_index_bp, adjusted_rand_index_bp, rand_index_seq, adjusted_rand_index_seq, \
           precision_avg_bp, precision_avg_seq, precision_weighted_bp, precision_weighted_seq, \
           recall_avg_bp, recall_avg_seq, recall_weighted_bp, recall_weighted_seq


def compute_metrics_tax(query, rank):
    print(query.label)
    gs_df = query.gold_standard_df[rank].reset_index()

    query_df = query.rank_to_df[rank].reset_index()

    query_w_length_df = pd.merge(query_df, gs_df.rename(columns={'TAXID': 'true_taxid'}).reset_index(), on='SEQUENCEID', sort=False)
    confusion_df = query_w_length_df.groupby(['TAXID', 'true_taxid'], sort=False).agg({'LENGTH': 'sum', 'SEQUENCEID': 'count'}).rename(columns={'LENGTH': 'tax_length', 'SEQUENCEID': 'tax_seq_counts'})
    rand_index_bp, adjusted_rand_index_bp = get_rand_index(confusion_df, 'TAXID', 'true_taxid', 'tax_length')
    rand_index_seq, adjusted_rand_index_seq = get_rand_index(confusion_df, 'TAXID', 'true_taxid', 'tax_seq_counts')

    query_w_length_df = query_w_length_df[['SEQUENCEID', 'TAXID', 'LENGTH']]

    percentage_of_assigned_bps = query_w_length_df['LENGTH'].sum() / gs_df['LENGTH'].sum()
    percentage_of_assigned_seqs = query_w_length_df.shape[0] / gs_df.shape[0]

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

    tp_length_sum = tp_fp_fn_df['tp_length'].sum()
    tp_seq_counts_sum = tp_fp_fn_df['tp_seq_counts'].sum()
    length_gs_sum = tp_fp_fn_df['length_gs'].sum()
    seq_counts_gs_sum = tp_fp_fn_df['seq_counts_gs'].sum()

    precision_avg_bp = tp_fp_fn_df['precision_bp'].mean()
    precision_avg_seq = tp_fp_fn_df['precision_seq'].mean()
    precision_weighted_bp = tp_length_sum / tp_fp_fn_df['total_length'].sum()
    precision_weighted_seq = tp_seq_counts_sum / tp_fp_fn_df['total_seq_counts'].sum()

    recall_avg_bp = tp_fp_fn_df['recall_bp'].mean()
    recall_avg_seq = tp_fp_fn_df['recall_seq'].mean()
    recall_weighted_bp = tp_length_sum / length_gs_sum
    recall_weighted_seq = tp_seq_counts_sum / seq_counts_gs_sum

    accuracy_bp = tp_length_sum / length_gs_sum
    accuracy_seq = tp_seq_counts_sum / seq_counts_gs_sum

    query.precision_df = tp_fp_fn_df
    query.recall_df = tp_fp_fn_df

    # f1_score_bp = 2 * precision_avg_bp * recall_avg_bp / (precision_avg_bp + recall_avg_bp)
    # f1_score_weighted_bp = 2 * precision_weighted_bp * recall_weighted_bp / (precision_weighted_bp + recall_weighted_bp)
    # print('f1_score_bp: {}'.format(f1_score_bp))

    return percentage_of_assigned_bps, percentage_of_assigned_seqs, accuracy_bp, accuracy_seq, \
           rand_index_bp, adjusted_rand_index_bp, rand_index_seq, adjusted_rand_index_seq, \
           precision_avg_bp, precision_avg_seq, precision_weighted_bp, precision_weighted_seq, \
           recall_avg_bp, recall_avg_seq, recall_weighted_bp, recall_weighted_seq


def evaluate(queries_list, sample_id):
    pd_bins_all = pd.DataFrame()
    df_summary = pd.DataFrame()

    for query in queries_list:

        if isinstance(query, new_classes.GenomeQueryNEW):

            percentage_of_assigned_bps, percentage_of_assigned_seqs, accuracy_bp, accuracy_seq, \
            rand_index_bp, adjusted_rand_index_bp, rand_index_seq, adjusted_rand_index_seq, \
            precision_avg_bp, precision_avg_seq, precision_weighted_bp, precision_weighted_seq, \
            recall_avg_bp, recall_avg_seq, recall_weighted_bp, recall_weighted_seq = compute_metrics(query)

        else:

            # for rank in query.rank_to_df:
            #     percentage_of_assigned_bps, percentage_of_assigned_seqs, accuracy_bp, accuracy_seq, \
            #     rand_index_bp, adjusted_rand_index_bp, rand_index_seq, adjusted_rand_index_seq, \
            #     precision_avg_bp, precision_avg_seq, precision_weighted_bp, precision_weighted_seq, \
            #     recall_avg_bp, recall_avg_seq, recall_weighted_bp, recall_weighted_seq = compute_metrics_tax(query, rank)
            percentage_of_assigned_bps, percentage_of_assigned_seqs, accuracy_bp, accuracy_seq, \
            rand_index_bp, adjusted_rand_index_bp, rand_index_seq, adjusted_rand_index_seq, \
            precision_avg_bp, precision_avg_seq, precision_weighted_bp, precision_weighted_seq, \
            recall_avg_bp, recall_avg_seq, recall_weighted_bp, recall_weighted_seq = compute_metrics_tax(query, 'species')


        precision_df = query.precision_df
        recall_df = query.recall_df

        df = pd.DataFrame(OrderedDict([(utils_labels.TOOL, query.label),
                                   (utils_labels.BINNING_TYPE, query.binning_type),
                                   (utils_labels.SAMPLE, sample_id),
                                   (utils_labels.RANK, 'NA'),
                                   (utils_labels.AVG_PRECISION_BP, [precision_avg_bp]),
                                   (utils_labels.AVG_PRECISION_BP_SEM, [precision_df['precision_bp'].sem()]),

                                   ('avg_precision_bp_var', [precision_df['precision_bp'].var()]),
                                   (utils_labels.AVG_RECALL_BP, [recall_avg_bp]),
                                   (utils_labels.AVG_RECALL_BP_SEM, [recall_df['recall_bp'].sem()]),
                                   ('avg_recall_bp_var', [recall_df['recall_bp'].var()]),
                                   (utils_labels.F1_SCORE_BP, [2 * precision_avg_bp * recall_avg_bp / (precision_avg_bp + recall_avg_bp)]),

                                   (utils_labels.AVG_PRECISION_SEQ, [precision_avg_seq]),
                                   (utils_labels.AVG_PRECISION_SEQ_SEM, [precision_df['precision_seq'].sem()]),
                                   (utils_labels.AVG_RECALL_SEQ, [recall_avg_seq]),
                                   (utils_labels.AVG_RECALL_SEQ_SEM, [recall_df['recall_seq'].sem()]),
                                   (utils_labels.F1_SCORE_SEQ, [2 * precision_avg_seq * recall_avg_seq / (precision_avg_seq + recall_avg_seq)]),

                                   (utils_labels.PRECISION_PER_BP, [precision_weighted_bp]),
                                   (utils_labels.PRECISION_PER_SEQ, [precision_weighted_seq]),
                                   (utils_labels.RECALL_PER_BP, [recall_weighted_bp]),
                                   (utils_labels.RECALL_PER_SEQ, [recall_weighted_seq]),
                                   (utils_labels.F1_SCORE_PER_BP, [2 * precision_weighted_bp * recall_weighted_bp / (precision_weighted_bp + recall_weighted_bp)]),
                                   (utils_labels.F1_SCORE_PER_SEQ, [2 * precision_weighted_seq * recall_weighted_seq / (precision_weighted_seq + recall_weighted_seq)]),

                                   (utils_labels.ACCURACY_PER_BP, [accuracy_bp]),
                                   (utils_labels.ACCURACY_PER_SEQ, [accuracy_seq]),

                                   (utils_labels.PERCENTAGE_ASSIGNED_BPS, [percentage_of_assigned_bps]),
                                   (utils_labels.PERCENTAGE_ASSIGNED_SEQS, [percentage_of_assigned_seqs]),
                                   (utils_labels.RI_BY_BP, [rand_index_bp]),
                                   (utils_labels.RI_BY_SEQ, [rand_index_seq]),
                                   (utils_labels.ARI_BY_BP, [adjusted_rand_index_bp]),
                                   (utils_labels.ARI_BY_SEQ, [adjusted_rand_index_seq]),

                                   (utils_labels.UNIFRAC_BP, [None]),
                                   (utils_labels.UNIFRAC_SEQ, [None]),

                                   (utils_labels.MISCLASSIFICATION_PER_BP, [1 - precision_weighted_bp]),
                                   (utils_labels.MISCLASSIFICATION_PER_SEQ, [1 - precision_weighted_seq])]))

        df_summary = pd.concat([df_summary, df], ignore_index=True, sort=True)

        precision_df[utils_labels.TOOL] = query.label
        pd_bins_all = pd.concat([pd_bins_all, precision_df.reset_index()], ignore_index=True, sort=False)

        print('percentage_of_assigned_bps: {}'.format(percentage_of_assigned_bps))
        print('percentage_of_assigned_seqs: {}'.format(percentage_of_assigned_seqs))
        print('accuracy_bp: {}'.format(accuracy_bp))
        print('accuracy_seq: {}'.format(accuracy_seq))

        print('rand_index_bp: {}'.format(rand_index_bp))
        print('adjusted_rand_index_bp: {}'.format(adjusted_rand_index_bp))
        print('rand_index_seq: {}'.format(rand_index_seq))
        print('adjusted_rand_index_seq: {}'.format(adjusted_rand_index_seq))

        print('precision_avg_bp: {}'.format(precision_avg_bp))
        print('precision_avg_seq: {}'.format(precision_avg_seq))
        print('precision_weighted_bp: {}'.format(precision_weighted_bp))
        print('precision_weighted_seq: {}'.format(precision_weighted_seq))

        print('recall_avg_bp: {}'.format(recall_avg_bp))
        print('recall_avg_seq: {}'.format(recall_avg_seq))
        print('recall_weighted_bp: {}'.format(recall_weighted_bp))
        print('recall_weighted_seq: {}'.format(recall_weighted_seq))

    pd_bins_all['sample_id'] = sample_id
    pd_bins_all['rank'] = 'NA'
    exit()

    return df_summary, pd_bins_all


def evaluate_samples_queries(sample_id_to_queries_list):
    pd_bins_all = pd.DataFrame()
    df_summary_all = pd.DataFrame()

    for sample_id in sample_id_to_queries_list:
        df_summary, pd_bins = evaluate(sample_id_to_queries_list[sample_id], sample_id)
        pd_bins_all = pd.concat([pd_bins_all, pd_bins], ignore_index=True)
        df_summary_all = pd.concat([df_summary_all, df_summary], ignore_index=True)
    return df_summary_all, pd_bins_all
