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
from concurrent.futures import ThreadPoolExecutor, wait
import pandas as pd
import os


def get_gs_data(query):
    return query.gold_standard_data


def compute_metrics(query, gs):
    query.compute_metrics(gs)


def evaluate_sample(pool, pool_io, queries_list):
    query1 = queries_list[0]
    gs_data = pool_io.submit(get_gs_data, query1).result()
    futures = [pool.submit(compute_metrics, query, gs_data) for query in queries_list]
    wait(futures)
    del gs_data


def evaluate_samples_queries(sample_id_to_g_queries_list, sample_id_to_t_queries_list, max_workers=None):
    if max_workers and max_workers <= 0:
        raise ValueError('Number of threads must be greater than 0')
    pool_main = ThreadPoolExecutor(max_workers)
    pool_child = ThreadPoolExecutor(max_workers)
    pool_io = ThreadPoolExecutor(min(10, os.cpu_count() or 1, max_workers)) if max_workers else ThreadPoolExecutor(min(10, os.cpu_count() or 1))

    for sample_id in sample_id_to_g_queries_list:
        if not sample_id_to_g_queries_list[sample_id]:
            continue
        pool_main.submit(evaluate_sample, pool_child, pool_io, sample_id_to_g_queries_list[sample_id])
    for sample_id in sample_id_to_t_queries_list:
        if not sample_id_to_g_queries_list[sample_id]:
            continue
        pool_main.submit(evaluate_sample, pool_child, pool_io, sample_id_to_t_queries_list[sample_id])
    pool_main.shutdown()
    pool_child.shutdown()
    pool_io.shutdown()

    pd_bins = pd.DataFrame()
    df_summary = pd.DataFrame()
    for sample_id in sample_id_to_g_queries_list:
        for query in sample_id_to_g_queries_list[sample_id]:
            if query.eval_success:
                df_summary = pd.concat([df_summary, query.get_metrics_df().dropna(axis=1, how='all')], ignore_index=True, sort=True)
                pd_bins = pd.concat([pd_bins, query.precision_df.reset_index()], ignore_index=True, sort=True)

    # Gold standard only has unfiltered metrics, so copy values to unfiltered columns
    for col in df_summary.columns:
        if col.endswith(utils_labels.UNFILTERED):
            df_summary.loc[df_summary[utils_labels.TOOL] == utils_labels.GS, col] = \
                df_summary.loc[df_summary[utils_labels.TOOL] == utils_labels.GS, col[:-len(utils_labels.UNFILTERED)]]

    return df_summary, pd_bins
