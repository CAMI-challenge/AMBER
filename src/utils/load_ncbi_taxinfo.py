#!/usr/bin/env python

RANKS = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain']
DICT_RANK_TO_INDEX = dict(zip(RANKS, list(range(len(RANKS)))))
RANKS_LOW2HIGH = list(reversed(RANKS))


def load_merged(names_file_path):
    tax_id_to_tax_id = {}
    with open(names_file_path) as read_handler:
        for line in read_handler:
            if len(line.strip()) == 0:
                continue
            line = line.split('|')
            line = list(map(str.strip, line))
            tax_id_to_tax_id[line[0]] = line[1]
    return tax_id_to_tax_id


def load_names(tax_id_to_rank, names_file_path):
    tax_id_to_name = {}
    with open(names_file_path) as read_handler:
        for line in read_handler:
            if len(line.strip()) == 0:
                continue
            line = line.split('|')
            line = list(map(str.strip, line))
            tax_id = line[0]

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
        elif tax_id_to_parent[tax_id] != '1' and tax_id_to_rank[tax_id_to_parent[tax_id]] not in RANKS:
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
            tax_id = line[0]
            tax_id_to_parent[tax_id] = line[1]
            tax_id_to_rank[tax_id] = line[2]

    for tax_id, rank in tax_id_to_rank.items():
        if tax_id_to_rank[tax_id] == "no rank" and check_parent_is_species(tax_id_to_parent, tax_id_to_rank, tax_id):
            tax_id_to_rank[tax_id] = "strain"

    return tax_id_to_parent, tax_id_to_rank


def get_id_path(tax_id, tax_id_to_parent, tax_id_to_rank):
    if tax_id not in tax_id_to_rank:
        # TODO report this in a log file
        return None

    while tax_id_to_rank[tax_id] not in RANKS:
        tax_id = tax_id_to_parent[tax_id]
        if tax_id == '1':
            return None

    index = DICT_RANK_TO_INDEX[tax_id_to_rank[tax_id]]
    path = [''] * (index + 1)
    path[index] = tax_id

    id = tax_id
    while id in tax_id_to_parent:
        id = tax_id_to_parent[id]
        if id == '1':
            break
        if tax_id_to_rank[id] not in RANKS:
            continue
        index = DICT_RANK_TO_INDEX[tax_id_to_rank[id]]
        path[index] = id
        if tax_id_to_rank[id] == "superkingdom":
            break
    return path
