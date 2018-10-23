#!/usr/bin/env python

import numpy as np
from abc import ABC, abstractmethod
from collections import defaultdict
from src.utils import load_ncbi_taxinfo
from src.utils import exclude_genomes
from src.utils import filter_tail


class Query(ABC):
    def __init__(self):
        self.__bins = []
        self.__bin_id_to_bin = {}
        self.__label = ""

    @property
    def bins(self):
        return self.__bins

    @property
    def label(self):
        return self.__label

    @bins.setter
    def bins(self, bins):
        self.__bins = bins

    @label.setter
    def label(self, label):
        self.__label = label

    def add_bin(self, bin):
        self.__bins.append(bin)
        self.__bin_id_to_bin[bin.id] = bin

    def get_bin_ids(self):
        return self.__bin_id_to_bin.keys()

    def get_sequence_ids(self):
        return set.union(*(bin.sequence_ids for bin in self.__bins))

    def get_bin_by_id(self, id):
        return self.__bin_id_to_bin[id]

    def get_bins_by_id(self, ids):
        return [self.get_bin_by_id(id) for id in ids]

    def get_all_mapping_ids(self):
        return [bin.mapping_id for bin in self.bins]


class GenomeQuery(Query):
    binning_type = 'genome'

    def __init__(self):
        super().__init__()
        self.__sequence_id_to_bin_id = {}
        self.__map_by_completeness = False

    @property
    def sequence_id_to_bin_id(self):
        return self.__sequence_id_to_bin_id

    @property
    def map_by_completeness(self):
        return self.__map_by_completeness

    @sequence_id_to_bin_id.setter
    def sequence_id_to_bin_id(self, sequence_id_bin_id):
        (sequence_id, bin_id) = sequence_id_bin_id
        self.__sequence_id_to_bin_id[sequence_id] = bin_id

    @map_by_completeness.setter
    def map_by_completeness(self, map_by_completeness):
        self.__map_by_completeness = map_by_completeness

    def compute_true_positives(self, gold_standard):
        for bin in self.bins:
            bin.compute_true_positives(gold_standard, self.__map_by_completeness)

    def get_bins_metrics(self, gold_standard):
        bins_metrics = [bin.get_metrics_dict(gold_standard) for bin in self.bins]
        mapped_ids = self.get_all_mapping_ids()
        for gs_bin in gold_standard.genome_query.bins:
            if gs_bin.id not in mapped_ids:
                bins_metrics.append({'id': None,
                                     'mapping_id': gs_bin.id,
                                     'purity': np.nan,
                                     'completeness': .0,
                                     'predicted_size': 0,
                                     'true_positives': 0,
                                     'real_size': gs_bin.length})

        if gold_standard.filter_tail_percentage:
            filter_tail.filter_tail(bins_metrics, gold_standard.filter_tail_percentage)
        if gold_standard.filter_genomes_file:
            bins_metrics = exclude_genomes.filter_data(bins_metrics, gold_standard.filter_genomes_file, gold_standard.filter_keyword)

        # sort bins by completeness
        return sorted(bins_metrics, key=lambda t: t['completeness'], reverse=True)


class TaxonomicQuery(Query):
    tax_id_to_rank = None
    binning_type = 'taxonomic'
    rank_to_overbinned_seqs = defaultdict(list)

    def __init__(self):
        super().__init__()
        self.__rank_to_sequence_id_to_bin_id = defaultdict(dict)

    @property
    def rank_to_sequence_id_to_bin_id(self):
        return self.__rank_to_sequence_id_to_bin_id

    @rank_to_sequence_id_to_bin_id.setter
    def rank_to_sequence_id_to_bin_id(self, rank_sequence_id_bin_id):
        (rank, sequence_id, bin_id) = rank_sequence_id_bin_id
        self.__rank_to_sequence_id_to_bin_id[rank][sequence_id] = bin_id

    def compute_true_positives(self, gold_standard):
        for bin in self.bins:
            bin.compute_true_positives(gold_standard)

    def get_bins_metrics(self, gold_standard):
        bins_metrics = []
        for bin in self.bins:
            bins_metrics.append(bin.get_metrics_dict(gold_standard))
        rank_to_index = dict(zip(load_ncbi_taxinfo.RANKS[::-1], list(range(len(load_ncbi_taxinfo.RANKS)))))
        # sort bins by rank and completeness
        return sorted(bins_metrics, key=lambda t: (rank_to_index[t['rank']], t['completeness']), reverse=True)


class Bin(ABC):
    sequence_id_to_length = {}

    def __init__(self, id):
        self.__id = id
        self.__sequence_ids = set()
        self.__length = 0
        self.__true_positives = 0
        self.__mapping_id = None
        self.__precision = .0
        self.__recall = .0
        self.__mapping_id_to_length = defaultdict(int)

    @property
    def id(self):
        return self.__id

    @property
    def sequence_ids(self):
        return self.__sequence_ids

    @property
    def length(self):
        return self.__length

    @property
    def true_positives(self):
        return self.__true_positives

    @property
    def mapping_id(self):
        return self.__mapping_id

    @property
    def precision(self):
        return self.__precision

    @property
    def recall(self):
        return self.__recall

    @property
    def mapping_id_to_length(self):
        return self.__mapping_id_to_length

    @id.setter
    def id(self, id):
        self.__id = id

    @sequence_ids.setter
    def sequence_ids(self, sequence_ids):
        self.__sequence_ids = sequence_ids

    @length.setter
    def length(self, length):
        self.__length = length

    @true_positives.setter
    def true_positives(self, true_positives):
        self.__true_positives = true_positives

    @mapping_id.setter
    def mapping_id(self, mapping_id):
        self.__mapping_id = mapping_id

    @precision.setter
    def precision(self, precision):
        self.__precision = precision

    @recall.setter
    def recall(self, recall):
        self.__recall = recall

    def add_sequence_id(self, sequence_id):
        if sequence_id not in self.__sequence_ids:
            self.__sequence_ids.add(sequence_id)
            self.__length += self.sequence_id_to_length[sequence_id]

    @abstractmethod
    def compute_confusion_matrix(self, gold_standard_query):
        pass

    @abstractmethod
    def get_metrics_dict(self):
        pass


class GenomeBin(Bin):
    def __init__(self, id):
        super().__init__(id)

    def compute_confusion_matrix(self, gold_standard):
        gold_standard_query = gold_standard.genome_query
        for sequence_id in self.sequence_ids:
            mapping_id = gold_standard_query.sequence_id_to_bin_id[sequence_id]
            self.mapping_id_to_length[mapping_id] += Bin.sequence_id_to_length[sequence_id]

    def compute_true_positives(self, gold_standard, map_by_completeness):
        if len(self.mapping_id_to_length) == 0:
            self.compute_confusion_matrix(gold_standard)
        gold_standard_query = gold_standard.genome_query
        if map_by_completeness:
            max_genome_percentage = .0
            best_gs_bin = gold_standard_query.bins[0]
            for gs_bin in gold_standard_query.bins:
                genome_percentage = self.mapping_id_to_length[gs_bin.id] / gs_bin.length
                if max_genome_percentage < genome_percentage:
                    max_genome_percentage = genome_percentage
                    best_gs_bin = gs_bin
                elif max_genome_percentage == genome_percentage and gs_bin.length > best_gs_bin.length:
                    best_gs_bin = gs_bin
            self.mapping_id = best_gs_bin.id
            self.true_positives = self.mapping_id_to_length[best_gs_bin.id]
        else:
            self.mapping_id = max(self.mapping_id_to_length, key=self.mapping_id_to_length.get)
            self.true_positives = self.mapping_id_to_length[self.mapping_id]

    def get_metrics_dict(self, gold_standard):
        gold_standard_query = gold_standard.genome_query
        return {'id': self.id,
                'mapping_id': self.mapping_id,
                'purity': self.precision,
                'completeness': self.recall,
                'predicted_size': self.length,
                'true_positives': self.true_positives,
                'real_size': gold_standard_query.get_bin_by_id(self.mapping_id).length}


class TaxonomicBin(Bin):
    def __init__(self, id):
        super().__init__(id)
        self.__rank = None

    @property
    def rank(self):
        return self.__rank

    @rank.setter
    def rank(self, rank):
        self.__rank = rank

    @property
    def mapping_id(self):
        return self.id

    def compute_confusion_matrix(self, gold_standard):
        gold_standard_query = gold_standard.taxonomic_query
        if self.rank not in gold_standard_query.rank_to_sequence_id_to_bin_id:
            return
        for sequence_id in self.sequence_ids:
            if sequence_id in gold_standard_query.rank_to_sequence_id_to_bin_id[self.rank]:
                mapping_id = gold_standard_query.rank_to_sequence_id_to_bin_id[self.rank][sequence_id]
                self.mapping_id_to_length[mapping_id] += Bin.sequence_id_to_length[sequence_id]

    def compute_true_positives(self, gold_standard):
        if len(self.mapping_id_to_length) == 0:
            self.compute_confusion_matrix(gold_standard)
        gold_standard_query = gold_standard.taxonomic_query
        if not self.id in gold_standard_query.get_bin_ids():
            return
        for gs_bin in gold_standard_query.bins:
            if self.id == gs_bin.id:
                commmon_seq_ids = self.sequence_ids & gs_bin.sequence_ids
                for sequence_id in commmon_seq_ids:
                    self.true_positives += Bin.sequence_id_to_length[sequence_id]

    def get_metrics_dict(self, gold_standard):
        gold_standard_query = gold_standard.taxonomic_query
        return {'id': self.id,
                'rank': gold_standard_query.tax_id_to_rank[self.id],
                'mapping_id': self.id,
                'purity': self.precision,
                'completeness': self.recall,
                'predicted_size': self.length,
                'true_positives': self.true_positives,
                'real_size': gold_standard_query.get_bin_by_id(
                    self.id).length if self.id in gold_standard_query.get_bin_ids() else np.nan}


class GoldStandard:
    def __init__(self, genome_query, taxonomic_query):
        self.__genome_query = genome_query
        self.__taxonomic_query = taxonomic_query
        self.__filter_tail_percentage = .0
        self.__filter_genomes_file = None
        self.__filter_keyword = None

    @property
    def genome_query(self):
        return self.__genome_query

    @property
    def taxonomic_query(self):
        return self.__taxonomic_query

    @property
    def filter_tail_percentage(self):
        return self.__filter_tail_percentage

    @property
    def filter_genomes_file(self):
        return self.__filter_genomes_file

    @property
    def filter_keyword(self):
        return self.__filter_keyword

    @genome_query.setter
    def genome_query(self, genome_query):
        self.__genome_query = genome_query

    @taxonomic_query.setter
    def taxonomic_query(self, taxonomic_query):
        self.__taxonomic_query = taxonomic_query

    @filter_tail_percentage.setter
    def filter_tail_percentage(self, filter_tail_percentage):
        self.__filter_tail_percentage = filter_tail_percentage

    @filter_genomes_file.setter
    def filter_genomes_file(self, filter_genomes_file):
        self.__filter_genomes_file = filter_genomes_file

    @filter_keyword.setter
    def filter_keyword(self, filter_keyword):
        self.__filter_keyword = filter_keyword
