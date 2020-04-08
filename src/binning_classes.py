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

import numpy as np
import logging
from abc import ABC, abstractmethod
from collections import defaultdict
from src.utils import load_ncbi_taxinfo
from src.utils import exclude_genomes
from src.utils import filter_tail
from src.utils import ProfilingTools as pf
from src import unifrac_distance as uf


class Query(ABC):
    def __init__(self):
        self.__sequence_id_to_length = None
        self.__bins = []
        self.__bin_id_to_bin = {}
        self.__label = ""
        self.__options = None
        self.__bins_metrics = None
        self.__gold_standard = None

    @property
    def sequence_id_to_length(self):
        return self.__sequence_id_to_length

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

    @sequence_id_to_length.setter
    def sequence_id_to_length(self, sequence_id_to_length):
        self.__sequence_id_to_length = sequence_id_to_length

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

    def compute_true_positives(self):
        for bin in self.bins:
            bin.compute_true_positives(self.gold_standard)

    def compute_precision_recall(self):
        for bin in self.__bins:
            bin.compute_precision_recall(self.__gold_standard)

    def compute_unifrac(self):
        return None, None


class GenomeQuery(Query):
    binning_type = 'genome'

    def __init__(self):
        super().__init__()
        self.__sequence_id_to_bin_id = {}
        self.__genome_to_recall_bp = {}
        self.__genome_to_recall_seq = {}
        self.__recall_bp = .0
        self.__recall_seq = .0
        self.__avg_recall_bp = .0
        self.__avg_recall_seq = .0

    @property
    def sequence_id_to_bin_id(self):
        return self.__sequence_id_to_bin_id

    @property
    def genome_to_recall_bp(self):
        return self.__genome_to_recall_bp

    @property
    def genome_to_recall_seq(self):
        return self.__genome_to_recall_seq

    @property
    def recall_bp(self):
        return self.__recall_bp

    @property
    def recall_seq(self):
        return self.__recall_seq

    @property
    def avg_recall_bp(self):
        return self.__avg_recall_bp

    @property
    def avg_recall_seq(self):
        return self.__avg_recall_seq

    @sequence_id_to_bin_id.setter
    def sequence_id_to_bin_id(self, sequence_id_bin_id):
        (sequence_id, bin_id) = sequence_id_bin_id

        # check if sequence is already in a bin
        if sequence_id in self.__sequence_id_to_bin_id:
            if isinstance(self.__sequence_id_to_bin_id[sequence_id], str):
                self.__sequence_id_to_bin_id[sequence_id] = {self.__sequence_id_to_bin_id[sequence_id], bin_id}
            elif isinstance(self.__sequence_id_to_bin_id[sequence_id], set):
                self.__sequence_id_to_bin_id[sequence_id].update(bin_id)
        else:
            self.__sequence_id_to_bin_id[sequence_id] = bin_id
        self.get_bin_by_id(bin_id).add_sequence_id(sequence_id, self.gold_standard.sequence_id_to_length[sequence_id])

    def compute_precision_recall(self):
        super().compute_precision_recall()
        self.compute_genome_to_recall()

    def compute_genome_to_recall(self):
        true_positives_recall_bp = 0
        true_positives_recall_seq = 0
        gs_bp = 0
        gs_seq = 0
        for i, gs_bin in enumerate(self.gold_standard.bins, start=1):
            if self.options.genome_to_unique_common:
                if gs_bin.id in self.options.genome_to_unique_common and (not self.options.filter_keyword or self.options.genome_to_unique_common[gs_bin.id] == self.options.filter_keyword):
                    continue

            max_bin_assigns_bp = 0
            max_bin_assigns_seq = 0
            best_bin = None
            for bin in self.bins:
                if bin.mapping_id_to_length[gs_bin.id] > max_bin_assigns_bp:
                    max_bin_assigns_bp = bin.mapping_id_to_length[gs_bin.id]
                    best_bin = bin
                if bin.mapping_id_to_num_seqs[gs_bin.id] > max_bin_assigns_seq:
                    max_bin_assigns_seq = bin.mapping_id_to_num_seqs[gs_bin.id]

            if max_bin_assigns_bp > 0:
                best_bin.represents_genome = True

            true_positives_recall_bp += max_bin_assigns_bp
            true_positives_recall_seq += max_bin_assigns_seq
            self.__genome_to_recall_bp[gs_bin.id] = max_bin_assigns_bp / gs_bin.length
            self.__genome_to_recall_seq[gs_bin.id] = max_bin_assigns_seq / gs_bin.num_seqs()

            gs_bp += gs_bin.length
            gs_seq += gs_bin.num_seqs()

            self.__avg_recall_bp = self.__avg_recall_bp + (self.__genome_to_recall_bp[gs_bin.id] - self.__avg_recall_bp) / i
            self.__avg_recall_seq = self.__avg_recall_seq + (self.__genome_to_recall_seq[gs_bin.id] - self.__avg_recall_seq) / i

        self.__recall_bp = true_positives_recall_bp / gs_bp
        self.__recall_seq = true_positives_recall_seq / gs_seq

    def get_all_mapping_ids(self):
        return [bin.most_complete_genome for bin in self.bins]

    def get_bins_metrics(self):
        if self.bins_metrics:
            return self.bins_metrics
        self.bins_metrics = [bin.get_metrics_dict(self.gold_standard) for bin in self.bins]
        mapped_ids = self.get_all_mapping_ids()
        for gs_bin in self.gold_standard.bins:
            if gs_bin.id not in mapped_ids:
                self.bins_metrics.append({'id': None,
                                          'rank': 'NA',
                                          'most_abundant_genome': gs_bin.most_abundant_genome,
                                          'most_complete_genome': gs_bin.most_complete_genome,
                                          'precision_bp': np.nan,
                                          'precision_seq': np.nan,
                                          'recall_bp': .0,
                                          'recall_seq': .0,
                                          'predicted_size': 0,
                                          'predicted_num_seqs': 0,
                                          'true_positive_bps': 0,
                                          'true_positive_seqs': 0,
                                          'true_size': gs_bin.length,
                                          'true_size_recall': gs_bin.length,
                                          'true_num_seqs': gs_bin.num_seqs(),
                                          'true_num_seqs_recall': gs_bin.num_seqs(),
                                          'represents_genome': True})

        if self.options.filter_tail_percentage and self != self.gold_standard:
            filter_tail.filter_tail(self.bins_metrics, self.options.filter_tail_percentage)
        if self.options.genome_to_unique_common:
            self.bins_metrics = exclude_genomes.filter_data(self.bins_metrics, self.options.genome_to_unique_common, self.options.filter_keyword)

        # sort bins by completeness
        self.bins_metrics = sorted(self.bins_metrics, key=lambda t: t['recall_bp'], reverse=True)
        return self.bins_metrics


class TaxonomicQuery(Query):
    tax_id_to_parent = None
    tax_id_to_rank = None
    tax_id_to_name = None
    tax_id_to_tax_id = None
    binning_type = 'taxonomic'

    def __init__(self):
        super().__init__()
        self.__rank_to_sequence_id_to_bin_id = defaultdict(dict)
        self.__rank_to_bins = defaultdict(list)
        self.__profile_bp = None
        self.__profile_seq = None

    @property
    def rank_to_sequence_id_to_bin_id(self):
        return self.__rank_to_sequence_id_to_bin_id

    @property
    def rank_to_bins(self):
        return self.__rank_to_bins

    @rank_to_sequence_id_to_bin_id.setter
    def rank_to_sequence_id_to_bin_id(self, rank_sequence_id_bin_id):
        (rank, sequence_id, bin_id) = rank_sequence_id_bin_id
        if sequence_id in self.__rank_to_sequence_id_to_bin_id[rank]:
            logging.getLogger('amber').warning("Sequence {} cannot be in multiple bins at rank {}.".format(sequence_id, rank))
        else:
            self.__rank_to_sequence_id_to_bin_id[rank][sequence_id] = bin_id
            self.get_bin_by_id(bin_id).add_sequence_id(sequence_id, self.gold_standard.sequence_id_to_length[sequence_id])

    def _create_profile(self, percentage_property):
        if not self.bins_metrics:
            self.get_bins_metrics()

        class Prediction:
            def __init__(self):
                pass
        profile = []
        for metrics in self.bins_metrics:
            prediction = Prediction()
            prediction.taxid = metrics['mapping_id']
            prediction.rank = metrics['rank']
            prediction.percentage = metrics[percentage_property]
            prediction.taxpath = '|'.join(load_ncbi_taxinfo.get_id_path(metrics['mapping_id'], TaxonomicQuery.tax_id_to_parent, TaxonomicQuery.tax_id_to_rank))
            prediction.taxpathsn = None
            profile.append(prediction)
        return profile

    @property
    def profile_bp(self):
        if self.__profile_bp:
            return self.__profile_bp
        self.__profile_bp = self._create_profile('true_positive_bps')
        return self.__profile_bp

    @property
    def profile_seq(self):
        if self.__profile_seq:
            return self.__profile_seq
        self.__profile_seq = self._create_profile('true_positive_seqs')
        return self.__profile_seq

    def add_bin(self, bin):
        self.rank_to_bins[bin.rank].append(bin)
        super().add_bin(bin)

    def get_bins_metrics(self):
        if self.bins_metrics:
            return self.bins_metrics
        self.bins_metrics = [bin.get_metrics_dict(self.gold_standard) for bin in self.bins]
        for rank in load_ncbi_taxinfo.RANKS:
            if rank in self.gold_standard.rank_to_bins:
                gs_rank_to_bin_ids = set([bin.id for bin in self.gold_standard.rank_to_bins[rank]])
            else:
                continue
            if rank in self.rank_to_bins:
                self_rank_to_bin_ids = set([bin.id for bin in self.rank_to_bins[rank]])
                ids_in_gs_but_not_in_self = gs_rank_to_bin_ids - self_rank_to_bin_ids
            else:
                ids_in_gs_but_not_in_self = gs_rank_to_bin_ids
            for bin_id in ids_in_gs_but_not_in_self:
                self.bins_metrics.append({'id': None,
                                          'name': TaxonomicQuery.tax_id_to_name[bin_id] if TaxonomicQuery.tax_id_to_name else np.nan,
                                          'rank': rank,
                                          'mapping_id': bin_id,
                                          'precision_bp': np.nan,
                                          'precision_seq': np.nan,
                                          'recall_bp': .0,
                                          'recall_seq': .0,
                                          'predicted_size': 0,
                                          'predicted_num_seqs': 0,
                                          'true_positive_bps': 0,
                                          'true_positive_seqs': 0,
                                          'true_size': self.gold_standard.get_bin_by_id(bin_id).length,
                                          'true_num_seqs': self.gold_standard.get_bin_by_id(bin_id).num_seqs()})

        if self.options.filter_tail_percentage and self != self.gold_standard:
            filter_tail.filter_tail(self.bins_metrics, self.options.filter_tail_percentage)

        rank_to_index = dict(zip(load_ncbi_taxinfo.RANKS[::-1], list(range(len(load_ncbi_taxinfo.RANKS)))))
        # sort bins by rank and completeness
        self.bins_metrics = sorted(self.bins_metrics, key=lambda t: (rank_to_index[t['rank']], t['recall_bp']), reverse=True)
        return self.bins_metrics

    def compute_unifrac(self):
        pf_profile_bp = pf.Profile(profile=self.profile_bp)
        gs_pf_profile_bp = pf.Profile(profile=self.gold_standard.profile_bp)
        pf_profile_seq = pf.Profile(profile=self.profile_seq)
        gs_pf_profile_seq = pf.Profile(profile=self.gold_standard.profile_seq)
        return uf.compute_unifrac(gs_pf_profile_bp, pf_profile_bp)[0], uf.compute_unifrac(gs_pf_profile_seq, pf_profile_seq)[0]


class Bin(ABC):
    def __init__(self, id):
        self.__id = id
        self.__sequence_ids = set()
        self.__length = 0
        self.__true_positive_bps = 0
        self.__true_positive_seqs = 0
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

    def add_sequence_id(self, sequence_id, length):
        if sequence_id not in self.__sequence_ids:
            self.__sequence_ids.add(sequence_id)
            self.__length += length

    @abstractmethod
    def compute_confusion_matrix(self, gold_standard):
        pass

    @abstractmethod
    def compute_precision_recall(self, gold_standard):
        pass

    @abstractmethod
    def get_metrics_dict(self):
        pass


class GenomeBin(Bin):
    def __init__(self, id):
        super().__init__(id)
        self.__most_abundant_genome = None
        self.__most_complete_genome = None
        self.__true_positive_bps_recall = 0
        self.__true_positive_seqs_recall = 0
        self.__represents_genome = False

    @property
    def most_abundant_genome(self):
        return self.__most_abundant_genome

    @property
    def most_complete_genome(self):
        return self.__most_complete_genome

    @property
    def true_positive_bps_recall(self):
        return self.__true_positive_bps_recall

    @property
    def true_positive_seqs_recall(self):
        return self.__true_positive_seqs_recall

    @property
    def represents_genome(self):
        return self.__represents_genome

    @most_abundant_genome.setter
    def most_abundant_genome(self, most_abundant_genome):
        self.__most_abundant_genome = most_abundant_genome

    @most_complete_genome.setter
    def most_complete_genome(self, most_complete_genome):
        self.__most_complete_genome = most_complete_genome

    @true_positive_bps_recall.setter
    def true_positive_bps_recall(self, true_positive_bps_recall):
        self.__true_positive_bps_recall = true_positive_bps_recall

    @true_positive_seqs_recall.setter
    def true_positive_seqs_recall(self, true_positive_seqs_recall):
        self.__true_positive_seqs_recall = true_positive_seqs_recall

    @represents_genome.setter
    def represents_genome(self, represents_genome):
        self.__represents_genome = represents_genome

    def compute_most_abundant_genome(self):
        self.__most_abundant_genome = max(self.mapping_id_to_length, key=self.mapping_id_to_length.get)

    def compute_most_complete_genome(self, gold_standard):
        max_genome_percentage = .0
        best_gs_bin = gold_standard.bins[0]
        for gs_bin in gold_standard.bins:
            if gs_bin.id in self.mapping_id_to_length:
                genome_percentage = self.mapping_id_to_length[gs_bin.id] / gs_bin.length
            else:
                genome_percentage = .0
            if max_genome_percentage < genome_percentage:
                max_genome_percentage = genome_percentage
                best_gs_bin = gs_bin
            elif max_genome_percentage == genome_percentage and gs_bin.length > best_gs_bin.length:
                best_gs_bin = gs_bin
        self.__most_complete_genome = best_gs_bin.id

    def compute_confusion_matrix(self, gold_standard):
        for sequence_id in self.sequence_ids:
            genome_id = gold_standard.sequence_id_to_bin_id[sequence_id]
            if isinstance(genome_id, str):
                self.mapping_id_to_length[genome_id] += gold_standard.sequence_id_to_length[sequence_id]
                self.mapping_id_to_num_seqs[genome_id] += 1
            elif isinstance(genome_id, set):
                for x in genome_id:
                    self.mapping_id_to_length[x] += gold_standard.sequence_id_to_length[sequence_id]
                    self.mapping_id_to_num_seqs[x] += 1

    def compute_true_positives(self, gold_standard):
        if len(self.mapping_id_to_length) == 0:
            self.compute_confusion_matrix(gold_standard)
            self.compute_most_abundant_genome()
            self.compute_most_complete_genome(gold_standard)
        self.true_positive_bps = self.mapping_id_to_length[self.most_abundant_genome]
        self.true_positive_seqs = self.mapping_id_to_num_seqs[self.most_abundant_genome]
        self.true_positive_bps_recall = self.mapping_id_to_length[self.most_complete_genome]
        self.true_positive_seqs_recall = self.mapping_id_to_num_seqs[self.most_complete_genome]

    def compute_precision_recall(self, gold_standard):
        self.precision_bp = self.true_positive_bps / self.length
        self.precision_seq = self.true_positive_seqs / self.num_seqs()
        self.recall_bp = self.true_positive_bps_recall / gold_standard.get_bin_by_id(self.most_complete_genome).length
        self.recall_seq = self.true_positive_seqs_recall / gold_standard.get_bin_by_id(self.most_complete_genome).num_seqs()

    def get_metrics_dict(self, gold_standard):
        return {'id': self.id,
                'rank': 'NA',
                'most_abundant_genome': self.most_abundant_genome,
                'most_complete_genome': self.most_complete_genome,
                'precision_bp': self.precision_bp,
                'precision_seq': self.precision_seq,
                'recall_bp': self.recall_bp,
                'recall_seq': self.recall_seq,
                'predicted_size': self.length,
                'predicted_num_seqs': self.num_seqs(),
                'true_positive_bps': self.true_positive_bps,
                'true_positive_seqs': self.true_positive_seqs,
                'true_size': gold_standard.get_bin_by_id(self.most_abundant_genome).length,
                'true_size_recall': gold_standard.get_bin_by_id(self.most_complete_genome).length,
                'true_num_seqs': gold_standard.get_bin_by_id(self.most_abundant_genome).num_seqs(),
                'true_num_seqs_recall': gold_standard.get_bin_by_id(self.most_complete_genome).num_seqs(),
                'represents_genome': self.represents_genome}


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
                self.mapping_id_to_length[mapping_id] += gold_standard.sequence_id_to_length[sequence_id]
                self.mapping_id_to_num_seqs[mapping_id] += 1

    def compute_precision_recall(self, gold_standard):
        self.precision_bp = self.true_positive_bps / self.length
        self.precision_seq = self.true_positive_seqs / self.num_seqs()
        if self.id in gold_standard.get_bin_ids():
            self.recall_bp = self.true_positive_bps / gold_standard.get_bin_by_id(self.id).length
            self.recall_seq = self.true_positive_seqs / gold_standard.get_bin_by_id(self.id).num_seqs()
        else:
            self.recall_bp = np.nan
            self.recall_seq = np.nan

    def compute_true_positives(self, gold_standard):
        if len(self.mapping_id_to_length) == 0:
            self.compute_confusion_matrix(gold_standard)
        if self.id not in gold_standard.get_bin_ids():
            return
        gs_bin = gold_standard.get_bin_by_id(self.id)
        commmon_seq_ids = self.sequence_ids & gs_bin.sequence_ids
        for sequence_id in commmon_seq_ids:
            self.true_positive_bps += gold_standard.sequence_id_to_length[sequence_id]
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
                'precision_bp': self.precision_bp,
                'precision_seq': self.precision_seq,
                'recall_bp': self.recall_bp,
                'recall_seq': self.recall_seq,
                'predicted_size': self.length,
                'predicted_num_seqs': self.num_seqs(),
                'true_positive_bps': self.true_positive_bps,
                'true_positive_seqs': self.true_positive_seqs,
                'true_size': true_size,
                'true_num_seqs': true_num_seqs}


class Options:
    def __init__(self, filter_tail_percentage, genome_to_unique_common, filter_keyword, min_length, rank_as_genome_binning):
        self.__filter_tail_percentage = float(filter_tail_percentage) if filter_tail_percentage else .0
        self.__genome_to_unique_common = genome_to_unique_common
        self.__filter_keyword = filter_keyword
        self.__min_length = int(min_length) if min_length else 0
        if rank_as_genome_binning and rank_as_genome_binning not in load_ncbi_taxinfo.RANKS:
            exit("Not a valid rank to assess taxonomic binning as genome binning (option --rank_as_genome_binning): " + rank_as_genome_binning)
        self.__rank_as_genome_binning = rank_as_genome_binning

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
