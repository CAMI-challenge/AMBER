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

from collections import OrderedDict

TOOL = "Tool"
RANK = "rank"
BINNING_TYPE = "binning type"
SAMPLE = "Sample"
GS = 'Gold standard'
UNFILTERED = ' (unfiltered)'

AVG_PRECISION_BP = 'precision_avg_bp'
AVG_PRECISION_BP_SEM = 'precision_avg_bp_sem'

AVG_PRECISION_SEQ = 'precision_avg_seq'
AVG_PRECISION_SEQ_SEM = 'precision_avg_seq_sem'

AVG_RECALL_BP = 'recall_avg_bp'
AVG_RECALL_BP_CAMI1 = 'recall_avg_bp_cami1'
AVG_RECALL_BP_SEM = 'recall_avg_bp_sem'
AVG_RECALL_BP_SEM_CAMI1 = 'recall_avg_bp_sem_cami1'

AVG_RECALL_SEQ = 'recall_avg_seq'
AVG_RECALL_SEQ_CAMI1 = 'recall_avg_seq_cami1'
AVG_RECALL_SEQ_SEM = 'recall_avg_seq_sem'
AVG_RECALL_SEQ_SEM_CAMI1 = 'recall_avg_seq_sem_cami1'

F1_SCORE_BP = 'f1_score_bp'
F1_SCORE_SEQ = 'f1_score_seq'

F1_SCORE_BP_CAMI1 = 'f1_score_bp_cami1'
F1_SCORE_SEQ_CAMI1 = 'f1_score_seq_cami1'

PRECISION_PER_BP = 'precision_weighted_bp'
PRECISION_PER_SEQ = 'precision_weighted_seq'

RECALL_PER_BP = 'recall_weighted_bp'
RECALL_PER_SEQ = 'recall_weighted_seq'

F1_SCORE_PER_BP = 'f1_score_per_bp'
F1_SCORE_PER_SEQ = 'f1_score_per_seq'

RI_BY_BP = 'rand_index_bp'
RI_BY_SEQ = 'rand_index_seq'
ARI_BY_BP = 'adjusted_rand_index_bp'
ARI_BY_SEQ = 'adjusted_rand_index_seq'

PERCENTAGE_ASSIGNED_BPS = 'percentage_of_assigned_bps'
PERCENTAGE_ASSIGNED_SEQS = 'percentage_of_assigned_seqs'

ACCURACY_PER_BP = 'accuracy_bp'
ACCURACY_PER_SEQ = 'accuracy_seq'
MISCLASSIFICATION_PER_BP = 'misclassification_bp'
MISCLASSIFICATION_PER_SEQ = 'misclassification_seq'

UNIFRAC_BP = 'unifrac_bp'
UNIFRAC_SEQ = 'unifrac_seq'

QUALITY_OF_BINS = "Quality of bins: all bins have the same weight"
QUALITY_OF_SAMPLE = "Quality for sample"

TOOLTIP_AVG_PRECISION_BP = "Average fraction of base pairs in each predicted genome bin coming from the most abundant genome in each bin. It ranges from 0 (worst) to 1 (best)."
TOOLTIP_AVG_PRECISION_BP_TAX = "Average fraction of true positive base pairs in each predicted bin. It ranges from 0 (worst) to 1 (best)."
TOOLTIP_AVG_PRECISION_BP_SEM = "Standard error of the mean of average purity (bp)."

TOOLTIP_AVG_PRECISION_SEQ = "Average fraction of sequences in each predicted genome bin coming from the most abundant genome (in base pairs) in each bin. It ranges from 0 (worst) to 1 (best)."
TOOLTIP_AVG_PRECISION_SEQ_TAX = "Average fraction of true positive sequences in each predicted bin. It ranges from 0 (worst) to 1 (best)."
TOOLTIP_AVG_PRECISION_SEQ_SEM = "Standard error of the mean of average purity (seq)."

TOOLTIP_AVG_RECALL_BP = "Average fraction of base pairs of each genome, each of which have been put together in one predicted genome bin. For each genome, the bin that includes most of the base pairs from that genome is evaluated. This metric ranges from 0 (worst) to 1 (best)."
TOOLTIP_AVG_RECALL_BP_TAX = "Average fraction of true positive base pairs in each predicted bin from the true bin size. It ranges from 0 (worst) to 1 (best)."

TOOLTIP_AVG_RECALL_BP_CAMI1 = "Average completeness per genome in base pairs, as computed in AMBER v1 (Meyer et al. GigaScience 2018) and the 1st CAMI Challenge (Sczyrba et al. Nature Methods 2017). It ranges from 0 (worst) to 1 (best)."
TOOLTIP_AVG_RECALL_BP_SEM = "Standard error of the mean of average completeness (bp)."
TOOLTIP_AVG_RECALL_SEQ = "Average fraction of sequences of each genome, each of which have been put together in one predicted genome bin. For each genome, the bin that includes most of the base pairs from that genome is evaluated. This metric ranges from 0 (worst) to 1 (best)."
TOOLTIP_AVG_RECALL_SEQ_TAX = "Average fraction of true positive sequences in each predicted bin from the true bin size in number of sequences. It ranges from 0 (worst) to 1 (best)."

TOOLTIP_AVG_RECALL_SEQ_CAMI1 = "Average completeness per genome in number of sequences, as computed in AMBER v1 (Meyer et al. GigaScience 2018) and the 1st CAMI Challenge (Sczyrba et al. Nature Methods 2017). It ranges from 0 (worst) to 1 (best)."
TOOLTIP_AVG_RECALL_SEQ_SEM = "Standard error of the mean of average completeness (seq)."

TOOLTIP_PRECISION_PER_BP = "Sum of base pairs coming from the most abundant genome in each predicted genome bin divided by the sum of base pairs in all predicted bins. It ranges from 0 (worst) to 1 (best)."
TOOLTIP_PRECISION_PER_BP_TAX = "Sum of true positive base pairs in the predicted bins divided by the sum of base pairs in all predicted bins. It ranges from 0 (worst) to 1 (best)."

TOOLTIP_PRECISION_PER_SEQ = "Sum of #sequences coming from the most abundant genome (in base pairs) in each predicted genome bin divided by the sum of #sequences in all predicted bins. It ranges from 0 (worst) to 1 (best)."
TOOLTIP_PRECISION_PER_SEQ_TAX = "Sum of true positive sequences in the predicted bins divided by the sum of #sequences in all predicted bins. It ranges from 0 (worst) to 1 (best)."

TOOLTIP_RECALL_PER_BP = "Sum of base pairs of each genome, each of which have been put together in one predicted genome bin, divided by the sum of base pairs of all underlying genomes. For each genome, the bin that includes most of the base pairs from that genome is considered. It ranges from 0 (worst) to 1 (best)."
TOOLTIP_RECALL_PER_BP_TAX = "Sum of true positive base pairs in the predicted bins divided by the sum of base pairs in all underlying taxa. It ranges from 0 (worst) to 1 (best)."

TOOLTIP_RECALL_PER_SEQ = "Sum of #sequences of each genome, each of which have been put together in one predicted genome bin, divided by the sum of #sequences of all underlying genomes. For each genome, the bin that includes most of the base pairs from that genome is considered. It ranges from 0 (worst) to 1 (best)."
TOOLTIP_RECALL_PER_SEQ_TAX = "Sum of true positive sequences in the predicted bins divided by the sum of #sequences in all underlying taxa. It ranges from 0 (worst) to 1 (best)."

TOOLTIP_F1_SCORE_BP = "Two times the product of average purity (bp) and average completeness (bp) divided by their sum. It ranges from 0 (worst) to 1 (best)."
TOOLTIP_F1_SCORE_SEQ = "Two times the product of average purity (seq) and average completeness (seq) divided by their sum. It ranges from 0 (worst) to 1 (best)."
TOOLTIP_F1_SCORE_PER_BP = "Two times the product of purity (bp) and completeness (bp) divided by their sum. It ranges from 0 (worst) to 1 (best)."
TOOLTIP_F1_SCORE_PER_SEQ = "Two times the product of purity (seq) and completeness (seq) divided by their sum. It ranges from 0 (worst) to 1 (best)."

TOOLTIP_RI_BY_BP = "Overall resolution in base pairs of the underlying ground truth genomes on the binned part of the sample. It ranges from 0 (worst) to 1 (best)."
TOOLTIP_RI_BY_BP_TAX = "Overall resolution in base pairs of the underlying ground truth taxa on the binned part of the sample. It ranges from 0 (worst) to 1 (best)."

TOOLTIP_RI_BY_SEQ = "Overall resolution in sequences of the underlying ground truth genomes on the binned part of the sample. It ranges from 0 (worst) to 1 (best)."
TOOLTIP_RI_BY_SEQ_TAX = "Overall resolution in sequences of the underlying ground truth taxa on the binned part of the sample. It ranges from 0 (worst) to 1 (best)."

TOOLTIP_ARI_BY_BP = "Rand index (bp) ajusted by the expected Rand index of a random clustering. It ranges from 0 (worst) to 1 (best)."
TOOLTIP_ARI_BY_SEQ = "Rand index (seq) ajusted by the expected Rand index of a random clustering. It ranges from 0 (worst) to 1 (best)."
TOOLTIP_PERCENTAGE_ASSIGNED_BPS = "Fraction of base pairs from the complete sample that have been assigned to predicted bins. It ranges from 0 (worst) to 1 (best)."
TOOLTIP_PERCENTAGE_ASSIGNED_SEQS = "Fraction of sequences from the complete sample that have been assigned to predicted bins. It ranges from 0 (worst) to 1 (best)."

TOOLTIP_ACCURACY_PER_BP = "Sum of base pairs coming from the most abundant genome in each predicted genome bin divided by the sum of base pairs in the complete sample. It ranges from 0 (worst) to 1 (best)."
TOOLTIP_ACCURACY_PER_BP_TAX = "Equivalent to completeness (bp), i.e. sum of true positive base pairs in the predicted bins divided by the sum of base pairs in all underlying taxa. It ranges from 0 (worst) to 1 (best)."

TOOLTIP_ACCURACY_PER_SEQ = "Sum of #sequences coming from the most abundant genome in each predicted genome bin divided by the sum of #sequences in the complete sample. It ranges from 0 (worst) to 1 (best)."
TOOLTIP_ACCURACY_PER_SEQ_TAX = "Equivalent to completeness (seq), i.e. sum of true positive sequences in the predicted bins divided by the sum of #sequences in all underlying taxa. It ranges from 0 (worst) to 1 (best)."

TOOLTIP_MISCLASSIFICATION_PER_BP = "1 - purity (bp)"
TOOLTIP_MISCLASSIFICATION_PER_SEQ = "1 - purity (seq)"

TOOLTIP_UNIFRAC_BP = "Tree-based measure of similarity between the true and predicted base pair assignments at all taxonomic ranks ranging from 0 (high similarity) to 16 (low similarity)."
TOOLTIP_UNIFRAC_SEQ = "Tree-based measure of similarity between the true and predicted base pair assignments at all taxonomic ranks ranging from 0 (high similarity) to 16 (low similarity)."

LABELS = {'precision_avg_bp': 'Average purity (bp)',
          'precision_avg_seq': 'Average purity (seq)',
          'precision_avg_bp_sem': 'Std error of av. purity (bp)',
          'precision_avg_seq_sem': 'Std error of av. purity (seq)',
          'recall_avg_bp': 'Average completeness (bp)',
          'recall_avg_seq': 'Average completeness (seq)',
          'recall_avg_bp_cami1': 'CAMI 1 average completeness (bp)',
          'recall_avg_seq_cami1': 'CAMI 1 average completeness (seq)',
          'recall_avg_bp_sem': 'Std error of av. completeness (bp)',
          'recall_avg_seq_sem': 'Std error of av. completeness (seq)',
          'recall_avg_bp_sem_cami1': 'CAMI 1 std error of av. completeness (bp)',
          'recall_avg_seq_sem_cami1': 'CAMI 1 std error of av. completeness (seq)',
          'f1_score_bp': 'F1 score (bp)',
          'f1_score_seq': 'F1 score (seq)',
          'f1_score_bp_cami1': 'CAMI 1 F1 score (bp)',
          'f1_score_seq_cami1': 'CAMI 1 F1 score (seq)',
          'f1_score_per_bp': 'F1 score for sample (bp)',
          'f1_score_per_seq': 'F1 score for sample (seq)',
          'accuracy_bp': 'Accuracy (bp)',
          'accuracy_seq': 'Accuracy (seq)',
          'misclassification_bp': 'Misclassification rate (bp)',
          'misclassification_seq': 'Misclassification rate (seq)',
          'precision_weighted_bp': 'Purity (bp)',
          'precision_weighted_seq': 'Purity (seq)',
          'recall_weighted_bp': 'Completeness (bp)',
          'recall_weighted_seq': 'Completeness (seq)',
          'rand_index_bp': 'Rand index (bp)',
          'rand_index_seq': 'Rand index (seq)',
          'adjusted_rand_index_bp': 'Adjusted Rand index (bp)',
          'adjusted_rand_index_seq': 'Adjusted Rand index (seq)',
          'percentage_of_assigned_bps': 'Percentage of binned bp',
          'percentage_of_assigned_seqs': 'Percentage of binned sequences',
          'unifrac_bp': 'UniFrac (bp)',
          'unifrac_seq': 'UniFrac (seq)'}


LABELS1 = [TOOL,
           ACCURACY_PER_BP,
           ACCURACY_PER_SEQ,
           ARI_BY_BP,
           ARI_BY_SEQ,
           AVG_RECALL_BP,
           AVG_RECALL_SEQ,
           AVG_PRECISION_BP,
           AVG_PRECISION_SEQ,
           F1_SCORE_BP_CAMI1,
           F1_SCORE_SEQ_CAMI1,
           AVG_RECALL_BP_CAMI1,
           AVG_RECALL_SEQ_CAMI1,
           AVG_RECALL_BP_SEM_CAMI1,
           AVG_RECALL_SEQ_SEM_CAMI1,
           RECALL_PER_BP,
           RECALL_PER_SEQ,
           F1_SCORE_BP,
           F1_SCORE_SEQ,
           F1_SCORE_PER_BP,
           F1_SCORE_PER_SEQ,
           MISCLASSIFICATION_PER_BP,
           MISCLASSIFICATION_PER_SEQ,
           PERCENTAGE_ASSIGNED_BPS,
           PERCENTAGE_ASSIGNED_SEQS,
           PRECISION_PER_BP,
           PRECISION_PER_SEQ,
           RI_BY_BP,
           RI_BY_SEQ,
           SAMPLE,
           AVG_RECALL_BP_SEM,
           AVG_RECALL_SEQ_SEM,
           AVG_PRECISION_BP_SEM,
           AVG_PRECISION_SEQ_SEM,
           UNIFRAC_BP,
           UNIFRAC_SEQ,
           BINNING_TYPE,
           RANK]


def get_genome_bins_columns():
    return OrderedDict([('BINID', 'Bin ID'),
                        ('genome_id', 'Most abundant genome'),
                        ('precision_bp', PRECISION_PER_BP),
                        ('recall_bp', RECALL_PER_BP),
                        ('total_length', 'Bin size (bp)'),
                        ('tp_length', 'True positives (bp)'),
                        ('length_gs', 'True size of most abundant genome (bp)'),
                        ('precision_seq', PRECISION_PER_SEQ),
                        ('recall_seq', RECALL_PER_SEQ),
                        ('total_seq_counts', 'Bin size (seq)'),
                        ('tp_seq_counts', 'True positives (seq)'),
                        ('seq_counts_gs', 'True size of most abundant genome (seq)')])


def get_tax_bins_columns():
    return OrderedDict([('TAXID', 'Taxon ID'),
                        ('name', 'Scientific name'),
                        ('rank', 'Taxonomic rank'),
                        ('precision_bp', PRECISION_PER_BP),
                        ('recall_bp', RECALL_PER_BP),
                        ('total_length', 'Bin size (bp)'),
                        ('tp_length', 'True positives (bp)'),
                        ('length_gs', 'True size (bp)'),
                        ('precision_seq', PRECISION_PER_SEQ),
                        ('recall_seq', RECALL_PER_SEQ),
                        ('total_seq_counts', 'Bin size (seq)'),
                        ('tp_seq_counts', 'True positives (seq)'),
                        ('seq_counts_gs', 'True size (seq)'),
                        ('filtered', 'Filtered')])
