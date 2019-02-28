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
        self.__options = None
        self.__bins_metrics = None
        self.__gold_standard = None

    @property
    def options(self):
        return self.__options

    @property
    def bins(self):
        return self.__bins

    @property
    def label(self):
        return self.__label

    @property
    def bins_metrics(self):
        return self.__bins_metrics

    @property
    def gold_standard(self):
        return self.__gold_standard

    @options.setter
    def options(self, options: 'Options'):
        self.__options = options

    @bins.setter
    def bins(self, bins):
        self.__bins = bins

    @label.setter
    def label(self, label):
        self.__label = label

    @bins_metrics.setter
    def bins_metrics(self, bins_metrics):
        self.__bins_metrics = bins_metrics

    @gold_standard.setter
    def gold_standard(self, gold_standard):
        self.__gold_standard = gold_standard

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

    def compute_precision_recall(self):
        for bin in self.__bins:
            bin.compute_precision_recall(self.__gold_standard)


class GenomeQuery(Query):
    binning_type = 'genome'

    def __init__(self):
        super().__init__()
        self.__sequence_id_to_bin_id = {}

    @property
    def sequence_id_to_bin_id(self):
        return self.__sequence_id_to_bin_id

    @sequence_id_to_bin_id.setter
    def sequence_id_to_bin_id(self, sequence_id_bin_id):
        (sequence_id, bin_id) = sequence_id_bin_id
        self.__sequence_id_to_bin_id[sequence_id] = bin_id

    def compute_true_positives(self):
        for bin in self.bins:
            bin.compute_true_positives(self.gold_standard, self.options.map_by_completeness)

    def get_bins_metrics(self):
        if self.bins_metrics:
            return self.bins_metrics
        self.bins_metrics = [bin.get_metrics_dict(self.gold_standard) for bin in self.bins]
        mapped_ids = self.get_all_mapping_ids()
        for gs_bin in self.gold_standard.bins:
            if gs_bin.id not in mapped_ids:
                self.bins_metrics.append({'id': None,
                                          'rank': 'NA',
                                          'mapping_id': gs_bin.id,
                                          'purity_bp': np.nan,
                                          'purity_seq': np.nan,
                                          'completeness_bp': .0,
                                          'completeness_seq': .0,
                                          'predicted_size': 0,
                                          'predicted_num_seqs': 0,
                                          'true_positive_bps': 0,
                                          'true_positive_seqs': 0,
                                          'true_size': gs_bin.length,
                                          'true_num_seqs': gs_bin.num_seqs()})

        if self.options.filter_tail_percentage:
            filter_tail.filter_tail(self.bins_metrics, self.options.filter_tail_percentage)
        if self.options.filter_genomes_file:
            self.bins_metrics = exclude_genomes.filter_data(self.bins_metrics, self.options.filter_genomes_file, self.options.filter_keyword)

        # sort bins by completeness
        self.bins_metrics = sorted(self.bins_metrics, key=lambda t: t['completeness_bp'], reverse=True)
        return self.bins_metrics


class TaxonomicQuery(Query):
    tax_id_to_parent = None
    tax_id_to_rank = None
    tax_id_to_name = None
    binning_type = 'taxonomic'

    def __init__(self):
        super().__init__()
        self.__rank_to_sequence_id_to_bin_id = defaultdict(dict)
        self.__rank_to_overbinned_seqs = defaultdict(list)
        self.__rank_to_bins = defaultdict(list)

    @property
    def rank_to_sequence_id_to_bin_id(self):
        return self.__rank_to_sequence_id_to_bin_id

    @property
    def rank_to_overbinned_seqs(self):
        return self.__rank_to_overbinned_seqs

    @property
    def rank_to_bins(self):
        return self.__rank_to_bins

    @rank_to_sequence_id_to_bin_id.setter
    def rank_to_sequence_id_to_bin_id(self, rank_sequence_id_bin_id):
        (rank, sequence_id, bin_id) = rank_sequence_id_bin_id
        self.__rank_to_sequence_id_to_bin_id[rank][sequence_id] = bin_id

    @rank_to_overbinned_seqs.setter
    def rank_to_overbinned_seqs(self, rank_overbinned_seqs):
        rank, overbinned_seqs = rank_overbinned_seqs
        self.__rank_to_overbinned_seqs[rank] = overbinned_seqs

    def add_bin(self, bin):
        self.rank_to_bins[bin.rank].append(bin)
        super().add_bin(bin)

    def append_overbinned_seq_id(self, rank, seq_id):
        self.__rank_to_overbinned_seqs[rank].append(seq_id)

    def compute_true_positives(self):
        for bin in self.bins:
            bin.compute_true_positives(self.gold_standard)

    def get_bins_metrics(self):
        if self.bins_metrics:
            return self.bins_metrics
        self.bins_metrics = [bin.get_metrics_dict(self.gold_standard) for bin in self.bins]
        for rank in load_ncbi_taxinfo.RANKS:
            if rank in self.gold_standard.rank_to_bins:
                gs_rank_ids = set([bin.id for bin in self.gold_standard.rank_to_bins[rank]])
            else:
                continue
            if rank in self.rank_to_bins:
                self_rank_ids = set([bin.id for bin in self.rank_to_bins[rank]])
                ids_in_gs_but_not_in_self = gs_rank_ids - self_rank_ids
            else:
                ids_in_gs_but_not_in_self = gs_rank_ids
            for bin_id in ids_in_gs_but_not_in_self:
                self.bins_metrics.append({'id': None,
                                          'name': TaxonomicQuery.tax_id_to_name[bin_id] if TaxonomicQuery.tax_id_to_name else np.nan,
                                          'rank': rank,
                                          'mapping_id': bin_id,
                                          'purity_bp': np.nan,
                                          'purity_seq': np.nan,
                                          'completeness_bp': .0,
                                          'completeness_seq': .0,
                                          'predicted_size': 0,
                                          'predicted_num_seqs': 0,
                                          'true_positive_bps': 0,
                                          'true_positive_seqs': 0,
                                          'true_size': self.gold_standard.get_bin_by_id(bin_id).length,
                                          'true_num_seqs': self.gold_standard.get_bin_by_id(bin_id).num_seqs()})
        rank_to_index = dict(zip(load_ncbi_taxinfo.RANKS[::-1], list(range(len(load_ncbi_taxinfo.RANKS)))))
        # sort bins by rank and completeness
        self.bins_metrics = sorted(self.bins_metrics, key=lambda t: (rank_to_index[t['rank']], t['completeness_bp']), reverse=True)
        return self.bins_metrics


class Bin(ABC):
    sequence_id_to_length = {}

    def __init__(self, id):
        self.__id = id
        self.__sequence_ids = set()
        self.__length = 0
        self.__true_positive_bps = 0
        self.__true_positive_seqs = 0
        self.__mapping_id = None
        self.__precision_bp = .0
        self.__precision_seq = .0
        self.__recall_bp = .0
        self.__recall_seq = .0
        self.__mapping_id_to_length = defaultdict(int)
        self.__mapping_id_to_num_seqs = defaultdict(int)

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
    def true_positive_bps(self):
        return self.__true_positive_bps

    @property
    def true_positive_seqs(self):
        return self.__true_positive_seqs

    @property
    def mapping_id(self):
        return self.__mapping_id

    @property
    def precision_bp(self):
        return self.__precision_bp

    @property
    def precision_seq(self):
        return self.__precision_seq

    @property
    def recall_bp(self):
        return self.__recall_bp

    @property
    def recall_seq(self):
        return self.__recall_seq

    @property
    def mapping_id_to_length(self):
        return self.__mapping_id_to_length

    @property
    def mapping_id_to_num_seqs(self):
        return self.__mapping_id_to_num_seqs

    @id.setter
    def id(self, id):
        self.__id = id

    @sequence_ids.setter
    def sequence_ids(self, sequence_ids):
        self.__sequence_ids = sequence_ids

    @length.setter
    def length(self, length):
        self.__length = length

    @true_positive_bps.setter
    def true_positive_bps(self, true_positive_bps):
        self.__true_positive_bps = true_positive_bps

    @true_positive_seqs.setter
    def true_positive_seqs(self, true_positive_seqs):
        self.__true_positive_seqs = true_positive_seqs

    @mapping_id.setter
    def mapping_id(self, mapping_id):
        self.__mapping_id = mapping_id

    @precision_bp.setter
    def precision_bp(self, precision_bp):
        self.__precision_bp = precision_bp

    @precision_seq.setter
    def precision_seq(self, precision_seq):
        self.__precision_seq = precision_seq

    @recall_bp.setter
    def recall_bp(self, recall_bp):
        self.__recall_bp = recall_bp

    @recall_seq.setter
    def recall_seq(self, recall_seq):
        self.__recall_seq = recall_seq

    def num_seqs(self):
        return len(self.__sequence_ids)

    def add_sequence_id(self, sequence_id):
        if sequence_id not in self.__sequence_ids:
            self.__sequence_ids.add(sequence_id)
            self.__length += self.sequence_id_to_length[sequence_id]

    @abstractmethod
    def compute_confusion_matrix(self, gold_standard):
        pass

    def compute_precision_recall(self, gold_standard):
        self.__precision_bp = self.__true_positive_bps / self.__length
        self.__precision_seq = self.__true_positive_seqs / self.num_seqs()
        if self.mapping_id in gold_standard.get_bin_ids():
            self.__recall_bp = self.__true_positive_bps / gold_standard.get_bin_by_id(self.mapping_id).length
            self.__recall_seq = self.__true_positive_seqs / gold_standard.get_bin_by_id(self.mapping_id).num_seqs()

    @abstractmethod
    def get_metrics_dict(self):
        pass


class GenomeBin(Bin):
    def __init__(self, id):
        super().__init__(id)

    def compute_confusion_matrix(self, gold_standard):
        for sequence_id in self.sequence_ids:
            mapping_id = gold_standard.sequence_id_to_bin_id[sequence_id]
            self.mapping_id_to_length[mapping_id] += Bin.sequence_id_to_length[sequence_id]
            self.mapping_id_to_num_seqs[mapping_id] += 1

    def compute_true_positives(self, gold_standard, map_by_completeness):
        if len(self.mapping_id_to_length) == 0:
            self.compute_confusion_matrix(gold_standard)
        if map_by_completeness:
            max_genome_percentage = .0
            best_gs_bin = gold_standard.bins[0]
            for gs_bin in gold_standard.bins:
                genome_percentage = self.mapping_id_to_length[gs_bin.id] / gs_bin.length
                if max_genome_percentage < genome_percentage:
                    max_genome_percentage = genome_percentage
                    best_gs_bin = gs_bin
                elif max_genome_percentage == genome_percentage and gs_bin.length > best_gs_bin.length:
                    best_gs_bin = gs_bin
            self.mapping_id = best_gs_bin.id
            self.true_positive_bps = self.mapping_id_to_length[best_gs_bin.id]
            self.true_positive_seqs = self.mapping_id_to_num_seqs[best_gs_bin.id]
        else:
            self.mapping_id = max(self.mapping_id_to_length, key=self.mapping_id_to_length.get)
            self.true_positive_bps = self.mapping_id_to_length[self.mapping_id]
            self.true_positive_seqs = self.mapping_id_to_num_seqs[self.mapping_id]

    def get_metrics_dict(self, gold_standard):
        return {'id': self.id,
                'rank': 'NA',
                'mapping_id': self.mapping_id,
                'purity_bp': self.precision_bp,
                'purity_seq': self.precision_seq,
                'completeness_bp': self.recall_bp,
                'completeness_seq': self.recall_seq,
                'predicted_size': self.length,
                'predicted_num_seqs': self.num_seqs(),
                'true_positive_bps': self.true_positive_bps,
                'true_positive_seqs': self.true_positive_seqs,
                'true_size': gold_standard.get_bin_by_id(self.mapping_id).length,
                'true_num_seqs': gold_standard.get_bin_by_id(self.mapping_id).num_seqs()}


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
        if self.rank not in gold_standard.rank_to_sequence_id_to_bin_id:
            return
        for sequence_id in self.sequence_ids:
            if sequence_id in gold_standard.rank_to_sequence_id_to_bin_id[self.rank]:
                mapping_id = gold_standard.rank_to_sequence_id_to_bin_id[self.rank][sequence_id]
                self.mapping_id_to_length[mapping_id] += Bin.sequence_id_to_length[sequence_id]
                self.mapping_id_to_num_seqs[mapping_id] += 1

    def compute_true_positives(self, gold_standard):
        if len(self.mapping_id_to_length) == 0:
            self.compute_confusion_matrix(gold_standard)
        if self.id not in gold_standard.get_bin_ids():
            return
        gs_bin = gold_standard.get_bin_by_id(self.id)
        commmon_seq_ids = self.sequence_ids & gs_bin.sequence_ids
        for sequence_id in commmon_seq_ids:
            self.true_positive_bps += Bin.sequence_id_to_length[sequence_id]
        self.true_positive_seqs = len(commmon_seq_ids)

    def get_metrics_dict(self, gold_standard):
        if self.id in gold_standard.get_bin_ids():
            true_size = gold_standard.get_bin_by_id(self.id).length
            true_num_seqs = gold_standard.get_bin_by_id(self.id).num_seqs()
        else:
            true_size = true_num_seqs = np.nan
        return {'id': self.id,
                'name': gold_standard.tax_id_to_name[self.id] if gold_standard.tax_id_to_name else np.nan,
                'rank': self.rank,
                'mapping_id': self.id,
                'purity_bp': self.precision_bp,
                'purity_seq': self.precision_seq,
                'completeness_bp': self.recall_bp,
                'completeness_seq': self.recall_seq,
                'predicted_size': self.length,
                'predicted_num_seqs': self.num_seqs(),
                'true_positive_bps': self.true_positive_bps,
                'true_positive_seqs': self.true_positive_seqs,
                'true_size': true_size,
                'true_num_seqs': true_num_seqs}


class Options:
    def __init__(self, filter_tail_percentage, filter_genomes_file, filter_keyword, map_by_completeness, min_length):
        self.__filter_tail_percentage = float(filter_tail_percentage) if filter_tail_percentage else .0
        self.__filter_genomes_file = filter_genomes_file
        self.__filter_keyword = filter_keyword
        self.__map_by_completeness = map_by_completeness
        self.__min_length = int(min_length) if min_length else 0

    @property
    def filter_tail_percentage(self):
        return self.__filter_tail_percentage

    @property
    def filter_genomes_file(self):
        return self.__filter_genomes_file

    @property
    def filter_keyword(self):
        return self.__filter_keyword

    @property
    def map_by_completeness(self):
        return self.__map_by_completeness

    @property
    def min_length(self):
        return self.__min_length

    @filter_tail_percentage.setter
    def filter_tail_percentage(self, filter_tail_percentage):
        self.__filter_tail_percentage = filter_tail_percentage

    @filter_genomes_file.setter
    def filter_genomes_file(self, filter_genomes_file):
        self.__filter_genomes_file = filter_genomes_file

    @filter_keyword.setter
    def filter_keyword(self, filter_keyword):
        self.__filter_keyword = filter_keyword

    @map_by_completeness.setter
    def map_by_completeness(self, map_by_completeness):
        self.__map_by_completeness = map_by_completeness

    @min_length.setter
    def min_length(self, min_length):
        self.__min_length = min_length
