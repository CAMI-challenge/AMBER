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

import collections
import os, sys, inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)
import pandas as pd


def transform_confusion_matrix2(query):
    gs_df = query.gold_standard_df.rename(columns={'bin_id': 'genome_id'})
    precision_df = query.precision_df
    confusion_df = query.confusion_df
    genomes_df = pd.DataFrame(gs_df['genome_id'].unique(), columns=['genome_id'])

    all_genomes_sorted_by_size = gs_df.groupby('genome_id', sort=False).agg({'seq_length': 'sum'}).sort_values(by='seq_length', axis=0, ascending=False).index.tolist()

    genomes_sorted_list = precision_df[['genome_id', 'genome_length']].sort_values(by='genome_length', ascending=False).groupby(['genome_id'], sort=False).first().index.tolist()
    unmapped_genomes = set(genomes_df['genome_id']) - set(genomes_sorted_list)
    unmapped_genomes = [genome_id for genome_id in all_genomes_sorted_by_size if genome_id in unmapped_genomes]
    genomes_sorted_list += unmapped_genomes

    unbinned_genomes_df = genomes_df[~genomes_df['genome_id'].isin(confusion_df.reset_index()['genome_id'])]
    heatmap_df = confusion_df.reset_index().append(unbinned_genomes_df, sort=False).fillna({'genome_length': 0, 'genome_seq_counts': 0})
    heatmap_df = heatmap_df.pivot(index='bin_id', columns='genome_id', values='genome_length')
    heatmap_df = heatmap_df.merge(precision_df['genome_length'], on='bin_id', sort=False).sort_values(by='genome_length', axis=0, ascending=False)

    gs_df = gs_df.set_index(['SEQUENCEID', 'genome_id'])
    query_df = query.df.set_index(['SEQUENCEID', 'genome_id'])
    difference_df = gs_df.index.difference(query_df.index)
    unassigned_sequences = gs_df.loc[difference_df].groupby('genome_id', sort=False).agg({'seq_length': 'sum'})

    heatmap_df = heatmap_df.append(unassigned_sequences.T, sort=True)
    heatmap_df = heatmap_df[genomes_sorted_list].fillna(0)

    # genomes_sorted_list = recall_df[['genome_id', 'genome_length']].sort_values(by='genome_length', ascending=False)['genome_id'].tolist()
    # genomes_sorted_list += list(set(genomes_df['genome_id']) - set(genomes_sorted_list))
    #
    # unbinned_genomes_df = genomes_df[~genomes_df['genome_id'].isin(confusion_df.reset_index()['genome_id'])]
    # heatmap_df = confusion_df.reset_index().append(unbinned_genomes_df, sort=False).fillna({'genome_length': 0, 'genome_seq_counts': 0})
    # heatmap_df = heatmap_df.pivot(index='bin_id', columns='genome_id', values='genome_length')
    # heatmap_df = heatmap_df.merge(precision_df['genome_length'], on='bin_id', sort=False).sort_values(by='genome_length', axis=0, ascending=False)
    # heatmap_df = heatmap_df[genomes_sorted_list].fillna(0)

    return heatmap_df


def transform_confusion_matrix(query):
    gold_standard = query.gold_standard

    bin_id_to_genome_id_to_length = {}
    for bin in query.bins:
        bin_id_to_genome_id_to_length[bin.id] = bin.mapping_id_to_length
    df_confusion = pd.DataFrame(bin_id_to_genome_id_to_length).T

    query_sequence_ids = set(query.sequence_id_to_bin_id.keys())
    gs_sequence_ids = set(gold_standard.sequence_id_to_bin_id.keys())
    genome_id_to_unassigned_bps = collections.Counter()
    for unassigned_seq_id in gs_sequence_ids - query_sequence_ids:
        genome_id = gold_standard.sequence_id_to_bin_id[unassigned_seq_id]
        genome_id_to_unassigned_bps[genome_id] += gold_standard.sequence_id_to_length[unassigned_seq_id]

    df_unassigned = pd.DataFrame.from_dict(genome_id_to_unassigned_bps, orient='index').rename(columns={0: 'unassigned'}).T
    table = df_confusion.append(df_unassigned, sort=False)
    table.fillna(0, inplace=True)
    # use log scale
    # table = table.applymap(np.log).fillna(0)

    bin_id_to_mapped_genome = {}
    for bin in query.bins:
        bin_id_to_mapped_genome[bin.id] = bin.most_abundant_genome

    # sort bins by the number of true positives (length of mapped genome within the bin)
    bin_id_to_mapped_genome_by_length = collections.OrderedDict(sorted(bin_id_to_mapped_genome.items(), key=lambda t: bin_id_to_genome_id_to_length[t[0]][t[1]], reverse=True))

    # sort genomes
    genome_order = []
    for bin_id in bin_id_to_mapped_genome_by_length:
        mapped_genome = bin_id_to_mapped_genome_by_length[bin_id]
        if mapped_genome not in genome_order:
            genome_order.append(mapped_genome)
    genome_order += list(set(table.columns.values.tolist()) - set(genome_order))
    for genome_id in genome_id_to_unassigned_bps.keys():
        if genome_id not in genome_order:
            genome_order.append(genome_id)

    table = table.loc[list(bin_id_to_mapped_genome_by_length.keys()) + ['unassigned'], genome_order]

    for genome_id in gold_standard.get_bin_ids():
        if genome_id not in table.columns.values.tolist():
            table[genome_id] = 0

    return table
