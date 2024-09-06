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

from cami_amber.utils import labels as utils_labels
import pandas as pd


def evaluate_sample(queries_list):
    query1 = queries_list[0]
    gs_data = query1.gold_standard_data
    for query in queries_list:
        query.compute_metrics(gs_data)


def evaluate_samples_queries(sample_id_to_g_queries_list, sample_id_to_t_queries_list):
    pd_bins = pd.DataFrame()
    df_summary = pd.DataFrame()

    def get_metrics(sample_id_to_queries_list):
        nonlocal pd_bins
        nonlocal df_summary

        for sample_id in sample_id_to_queries_list:
            evaluate_sample(sample_id_to_queries_list[sample_id])

        for sample_id in sample_id_to_queries_list:
            for query in sample_id_to_queries_list[sample_id]:
                if query.eval_success:
                    df_summary = pd.concat([df_summary, query.get_metrics_df().dropna(axis=1, how='all')], ignore_index=True, sort=True)
                    pd_bins = pd.concat([pd_bins, query.precision_df.reset_index()], ignore_index=True, sort=True)

    get_metrics(sample_id_to_g_queries_list)
    get_metrics(sample_id_to_t_queries_list)

    # Gold standard only has unfiltered metrics, so copy values to unfiltered columns
    for col in df_summary.columns:
        if col.endswith(utils_labels.UNFILTERED):
            df_summary.loc[df_summary[utils_labels.TOOL] == utils_labels.GS, col] = \
                df_summary.loc[df_summary[utils_labels.TOOL] == utils_labels.GS, col[:-len(utils_labels.UNFILTERED)]]

    return df_summary, pd_bins
