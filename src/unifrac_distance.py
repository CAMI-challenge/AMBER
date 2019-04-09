# This is an example of how to use the script
from src.utils import ProfilingTools as PF
from src.utils import EMDUnifrac as EMDU
import numpy as np
import argparse
import copy
import sys
import os


def compute_unifrac(pf1, pf2):
    P1 = copy.deepcopy(pf1)
    P2 = copy.deepcopy(pf2)
    (Tint, lint, nodes_in_order, nodes_to_index, P, Q) = P1.make_unifrac_input_and_normalize(P2)
    (weighted, _) = EMDU.EMDUnifrac_weighted(Tint, lint, nodes_in_order, P, Q)
    (unweighted, _) = EMDU.EMDUnifrac_unweighted(Tint, lint, nodes_in_order, P, Q)
    return weighted, unweighted


def weighted_unifrac(pf1, pf2):
    P1 = copy.deepcopy(pf1)
    P2 = copy.deepcopy(pf2)
    (Tint, lint, nodes_in_order, nodes_to_index, P, Q) = P1.make_unifrac_input_and_normalize(P2)
    (val, _) = EMDU.EMDUnifrac_weighted(Tint, lint, nodes_in_order, P, Q)
    return val


def unweighted_unifrac(pf1, pf2):
    P1 = copy.deepcopy(pf1)
    P2 = copy.deepcopy(pf2)
    (Tint, lint, nodes_in_order, nodes_to_index, P, Q) = P1.make_unifrac_input_and_normalize(P2)
    (val, _) = EMDU.EMDUnifrac_unweighted(Tint, lint, nodes_in_order, P, Q)
    return val


def print_list_unifrac(unifrac_list, labels, stream=sys.stdout):
    labels_iterator = iter(labels)
    stream.write("tool\tweighted\tunweighted\n")
    for unifrac in unifrac_list:
        stream.write("{}\t".format(next(labels_iterator)))
        stream.write("\t".join([format(unifrac[0], '.3f'), format(unifrac[1], '.3f')]))
        stream.write("\n")


def read_params(args):
	parser = argparse.ArgumentParser(description='')
	arg = parser.add_argument
	arg('--input', metavar='files_file', type=str, required=True,
					default=None, help="File of CAMI profile files to use")
	arg('--output', metavar='output_file', required=True, default=None, type=str,
					help="Output file (you should have this end in .csv as it is a matrix)")
	arg('--threshold', metavar='threshold', required=True, default=None, type=float,
					help="Value to threshold profiles to before computing EMDUnifrac. "
					"NOTE THIS VALUE IS IN PERCENTAGES so if you want 1% use 1")
	return vars(parser.parse_args())


if __name__ == '__main__':
	par = read_params(sys.argv)
	files_file = par['input']
	output_file = par['output']
	threshold = par['threshold']

	# Get all the profile file names
	files = []
	fid = open(files_file, 'r')
	for line in fid.readlines():
		files.append(line.strip())
	fid.close()

	# Import all the profiles
	profiles = []
	for file_name in files:
		profiles.append(PF.Profile(input_file_name=file_name))

	# Threshold all the profiles
	for profile in profiles:
		profile.threshold(threshold=threshold)

	# Compute EMDUnifrac
	D = np.zeros((len(profiles), len(profiles)))
	for i in range(len(profiles)):
		for j in range(i+1, len(profiles)):
			(wu, uu) = compute_unifrac(profiles[i], profiles[j])
			D[i, j] = wu
			D[j, i] = uu

	# Save results in tsv
	np.savetxt(output_file, D, delimiter='\t', newline='\n')
