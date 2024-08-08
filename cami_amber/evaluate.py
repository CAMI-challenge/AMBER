# Copyright 2023 Department of Computational Biology for Infection Research - Helmholtz Centre for Infection Research
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

from cami_amber.utils import labels as utils_labels
import pandas as pd


def evaluate_each(queries_list, sample_id, gs):
    pd_bins_all = pd.DataFrame()
    df_summary = pd.DataFrame()

    for query in queries_list:
        if not query.compute_metrics(gs):
            continue

        query_metrics_df = query.get_metrics_df()
        query_metrics_df[utils_labels.SAMPLE] = sample_id

        df_summary = pd.concat([df_summary, query_metrics_df.dropna(axis=1, how='all')], ignore_index=True, sort=True)

        query.precision_df[utils_labels.TOOL] = query.label
        pd_bins_all = pd.concat([pd_bins_all, query.precision_df.reset_index()], ignore_index=True, sort=True)

    pd_bins_all['sample_id'] = sample_id

    return df_summary, pd_bins_all


def evaluate_samples_queries(sample_id_to_g_queries_list, sample_id_to_t_queries_list):
    pd_bins_all = pd.DataFrame()
    df_summary_all = pd.DataFrame()

    for sample_id in sample_id_to_g_queries_list:
        if not sample_id_to_g_queries_list[sample_id]:
            continue
        query1 = sample_id_to_g_queries_list[sample_id][0]
        gs_df = query1.gold_standard.df
        gs_df = gs_df[['SEQUENCEID', 'BINID', 'LENGTH']].rename(columns={'LENGTH': 'seq_length', 'BINID': 'genome_id'})
        df_summary, pd_bins = evaluate_each(sample_id_to_g_queries_list[sample_id], sample_id, gs_df)
        del gs_df
        pd_bins_all = pd.concat([pd_bins_all, pd_bins], ignore_index=True)
        df_summary_all = pd.concat([df_summary_all, df_summary], ignore_index=True)

    for sample_id in sample_id_to_t_queries_list:
        if not sample_id_to_t_queries_list[sample_id]:
            continue
        query1 = sample_id_to_t_queries_list[sample_id][0]
        gs_rank_to_df = query1.gold_standard.rank_to_df
        df_summary, pd_bins = evaluate_each(sample_id_to_t_queries_list[sample_id], sample_id, gs_rank_to_df)
        del gs_rank_to_df
        pd_bins_all = pd.concat([pd_bins_all, pd_bins], ignore_index=True)
        df_summary_all = pd.concat([df_summary_all, df_summary], ignore_index=True)

    # Gold standard only has unfiltered metrics, so copy values to unfiltered columns
    for col in df_summary_all.columns:
        if col.endswith(utils_labels.UNFILTERED):
            df_summary_all.loc[df_summary_all[utils_labels.TOOL] == utils_labels.GS, col] = \
                df_summary_all.loc[df_summary_all[utils_labels.TOOL] == utils_labels.GS, col[:-len(utils_labels.UNFILTERED)]]

    return df_summary_all, pd_bins_all