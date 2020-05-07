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

TOOL = "Tool"
RANK = "rank"
BINNING_TYPE = "binning type"
SAMPLE = "Sample"
GS = 'Gold standard'

AVG_PRECISION_BP = "Average purity (bp)"
AVG_PRECISION_BP_SEM = "Std error of average purity (bp)"

AVG_PRECISION_SEQ = "Average purity (seq)"
AVG_PRECISION_SEQ_SEM = "Std error of average purity (seq)"

AVG_RECALL_BP = "Average completeness (bp) "
AVG_RECALL_BP_CAMI1 = "CAMI 1 average completeness (bp)"
AVG_RECALL_BP_SEM = "Std error of average completeness (bp)"

AVG_RECALL_SEQ = "Average completeness (seq) "
AVG_RECALL_SEQ_CAMI1 = "CAMI 1 average completeness (seq)"
AVG_RECALL_SEQ_SEM = "Std error of average completeness (seq)"

F1_SCORE_BP = "F1 score (bp)"
F1_SCORE_SEQ = "F1 score (seq)"

PRECISION_PER_BP = "Purity (bp)"
PRECISION_PER_SEQ = "Purity (seq)"

RECALL_PER_BP = "Completeness (bp)"
RECALL_PER_SEQ = "Completeness (seq)"

F1_SCORE_PER_BP = "F1 score for sample (bp)"
F1_SCORE_PER_SEQ = "F1 score for sample (seq)"

RI_BY_BP = "Rand index (bp)"
RI_BY_SEQ = "Rand index (seq)"
ARI_BY_BP = "Adjusted Rand index (bp)"
ARI_BY_SEQ = "Adjusted Rand index (seq)"

PERCENTAGE_ASSIGNED_BPS = "Percentage of binned bp"
PERCENTAGE_ASSIGNED_SEQS = "Percentage of binned sequences"

ACCURACY_PER_BP = "Accuracy (bp)"
ACCURACY_PER_SEQ = "Accuracy (seq)"
MISCLASSIFICATION_PER_BP = "Misclassification rate (bp)"
MISCLASSIFICATION_PER_SEQ = "Misclassification rate (seq)"

UNIFRAC_BP = "UniFrac (bp)"
UNIFRAC_SEQ = "UniFrac (seq)"

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

TOOLTIP_PRECISION_PER_SEQ = "Sum of sequences coming from the most abundant genome (in base pairs) in each predicted genome bin divided by the sum of sequences in all predicted bins. It ranges from 0 (worst) to 1 (best)."
TOOLTIP_PRECISION_PER_SEQ_TAX = "Sum of true positive sequences in the predicted bins divided by the sum of sequences in all predicted bins. It ranges from 0 (worst) to 1 (best)."

TOOLTIP_RECALL_PER_BP = "Sum of base pairs of each genome, each of which have been put together in one predicted genome bin, divided by the sum of base pairs of all underlying genomes. For each genome, the bin that includes most of the base pairs from that genome is considered. It ranges from 0 (worst) to 1 (best)."
TOOLTIP_RECALL_PER_BP_TAX = "Sum of true positive base pairs in the predicted bins divided by the sum of base pairs in all underlying taxa. It ranges from 0 (worst) to 1 (best)."

TOOLTIP_RECALL_PER_SEQ = "Sum of sequences of each genome, each of which have been put together in one predicted genome bin, divided by the sum of sequences of all underlying genomes. For each genome, the bin that includes most of the base pairs from that genome is considered. It ranges from 0 (worst) to 1 (best)."
TOOLTIP_RECALL_PER_SEQ_TAX = "Sum of true positive sequences in the predicted bins divided by the sum of sequences in all underlying taxa. It ranges from 0 (worst) to 1 (best)."

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

TOOLTIP_ACCURACY_PER_BP = "Sum of base pairs coming from the most abundant genome in each predicted genome bin divided by the sum of base pairs in the complete data set. It ranges from 0 (worst) to 1 (best)."
TOOLTIP_ACCURACY_PER_BP_TAX = "Equivalent to completeness (bp), i.e. sum of true positive base pairs in the predicted bins divided by the sum of sequences in all underlying taxa. It ranges from 0 (worst) to 1 (best)."

TOOLTIP_ACCURACY_PER_SEQ = "Sum of sequences coming from the most abundant genome in each predicted genome bin divided by the sum of sequences in the complete data set. It ranges from 0 (worst) to 1 (best)."
TOOLTIP_ACCURACY_PER_SEQ_TAX = "Equivalent to completeness (seq), i.e. sum of true positive sequences in the predicted bins divided by the sum of sequences in all underlying taxa. It ranges from 0 (worst) to 1 (best)."

TOOLTIP_MISCLASSIFICATION_PER_BP = "1 - purity (bp)"
TOOLTIP_MISCLASSIFICATION_PER_SEQ = "1 - purity (seq)"

TOOLTIP_UNIFRAC_BP = "Tree-based measure of similarity between the true and predicted base pair assignments at all taxonomic ranks ranging from 0 (high similarity) to 16 (low similarity)."
TOOLTIP_UNIFRAC_SEQ = "Tree-based measure of similarity between the true and predicted base pair assignments at all taxonomic ranks ranging from 0 (high similarity) to 16 (low similarity)."
