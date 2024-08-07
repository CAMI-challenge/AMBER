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

import logging
import os
import numpy as np
import pandas as pd

RANKS = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain']
RANK_TO_INDEX = dict(zip(RANKS, list(range(len(RANKS)))))
INDEX_TO_RANK = dict(zip(list(range(len(RANKS))), RANKS))
RANKS_LOW2HIGH = list(reversed(RANKS))


def load_merged(names_file_path):
    tax_id_to_tax_id = {}
    with open(names_file_path) as read_handler:
        for line in read_handler:
            if len(line.strip()) == 0:
                continue
            line = line.split('|')
            line = list(map(str.strip, line))
            tax_id_to_tax_id[int(line[0])] = int(line[1])
    return tax_id_to_tax_id


def load_names(tax_id_to_rank, names_file_path):
    tax_id_to_name = {}
    with open(names_file_path) as read_handler:
        for line in read_handler:
            if len(line.strip()) == 0:
                continue
            line = line.split('|')
            line = list(map(str.strip, line))
            tax_id = int(line[0])

            if line[3] == "scientific name":
                tax_id_to_name[tax_id] = line[1]
            else:
                continue

            if tax_id_to_rank[tax_id] == "species" or tax_id_to_rank[tax_id] == "strain":
                names = tax_id_to_name[tax_id].split(" ")
                if len(names) > 2:
                    if tax_id_to_rank[tax_id] == "strain":
                        tax_id_to_name[tax_id] = "{} {} strain".format(names[0], names[1])
                    else:
                        tax_id_to_name[tax_id] = "{} {}".format(names[0], names[1])
    tax_id_to_name[""] = ""
    return tax_id_to_name


def check_parent_is_species(tax_id_to_parent, tax_id_to_rank, tax_id):
    if tax_id_to_parent[tax_id] in tax_id_to_parent:
        if tax_id_to_rank[tax_id_to_parent[tax_id]] == "species":
            return True
        elif tax_id_to_parent[tax_id] != 1 and tax_id_to_rank[tax_id_to_parent[tax_id]] not in RANKS:
            return check_parent_is_species(tax_id_to_parent, tax_id_to_rank, tax_id_to_parent[tax_id])
    return False


def load_tax_info(ncbi_nodes_file):
    tax_id_to_parent = {}
    tax_id_to_rank = {}
    with open(ncbi_nodes_file) as read_handler:
        for line in read_handler:
            if len(line.strip()) == 0:
                continue
            line = line.split('|')
            line = list(map(str.strip, line))
            tax_id = int(line[0])
            tax_id_to_parent[tax_id] = int(line[1])
            tax_id_to_rank[tax_id] = line[2]

    for tax_id, rank in tax_id_to_rank.items():
        if tax_id_to_rank[tax_id] == "no rank" and check_parent_is_species(tax_id_to_parent, tax_id_to_rank, tax_id):
            tax_id_to_rank[tax_id] = "strain"

    return tax_id_to_parent, tax_id_to_rank


def get_id_path(tax_id, tax_id_to_parent, tax_id_to_rank, tax_id_to_tax_id):
    if tax_id not in tax_id_to_rank:
        if tax_id_to_tax_id and tax_id in tax_id_to_tax_id:
            tax_id = tax_id_to_tax_id[tax_id]
        else:
            logging.getLogger('amber').warning("Invalid NCBI taxonomic ID: {}".format(tax_id))
            return [np.nan] * 8

    while tax_id_to_rank[tax_id] not in RANKS:
        tax_id = tax_id_to_parent[tax_id]
        if tax_id == 1:
            return [np.nan] * 8

    index = RANK_TO_INDEX[tax_id_to_rank[tax_id]]
    path = [np.nan] * 8
    path[index] = tax_id

    id = tax_id
    try:
        while id in tax_id_to_parent:
            id = tax_id_to_parent[id]
            if id == 1:
                break
            if tax_id_to_rank[id] not in RANKS:
                continue
            index = RANK_TO_INDEX[tax_id_to_rank[id]]
            path[index] = id
            if tax_id_to_rank[id] == 'superkingdom':
                break
    except KeyError:
        logging.getLogger('amber').warning('Could not get rank for taxonomic ID: {}'.format(tax_id))
    return path


def preprocess_ncbi_tax(ncbi_dir):
    def get_path(row):
        return get_id_path(row['TAXID'], tax_id_to_parent, tax_id_to_rank, tax_id_to_tax_id)

    merged = pd.read_csv(os.path.join(ncbi_dir, 'merged.dmp'), sep='|', engine='python', header=None, index_col=False, names=['TAXID'], usecols=[0])
    nodes = pd.read_csv(os.path.join(ncbi_dir, 'nodes.dmp'), sep='|', engine='python', header=None, names=['TAXID'], usecols=[0])
    nodes = pd.concat([nodes, merged], ignore_index=True)

    tax_id_to_parent, tax_id_to_rank = load_tax_info(os.path.join(ncbi_dir, 'nodes.dmp'))
    tax_id_to_tax_id = load_merged(os.path.join(ncbi_dir, 'merged.dmp'))

    nodes = pd.concat([nodes, nodes.apply(lambda row: get_path(row), axis=1, result_type='expand')], axis='columns')
    for i in range(8):
        nodes[i] = nodes[i].astype(pd.UInt32Dtype())
    nodes.rename(columns=INDEX_TO_RANK, inplace=True)

    tax_id_to_name = load_names(tax_id_to_rank, os.path.join(ncbi_dir, 'names.dmp'))
    for k, v in tax_id_to_tax_id.items():
        tax_id_to_name[k] = tax_id_to_name[tax_id_to_tax_id[k]]

    nodes['name'] = nodes['TAXID'].map(tax_id_to_name)

    for k, v in tax_id_to_tax_id.items():
        tax_id_to_rank[k] = tax_id_to_rank[tax_id_to_tax_id[k]]
    nodes['rank'] = nodes.apply(lambda x: tax_id_to_rank[x['TAXID']], axis=1)

    if os.access(ncbi_dir, os.W_OK):
        nodes.to_feather(os.path.join(ncbi_dir, 'nodes.amber.ft'))
    else:
        logging.getLogger('amber').warning('You have no permission to write processed NCBI taxonomy to {}'.format(ncbi_dir))

    return nodes.set_index('TAXID')
