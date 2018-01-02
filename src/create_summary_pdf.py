#!/usr/bin/env python

import argparse
import os, sys, inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)
from src.utils import argparse_parents
from src.utils import load_data

PLOTS = ['avg_purity_completeness',
         'ari_vs_assigned_bps',
         'avg_purity_completeness_per_bp',
         'completeness_contamination_table',
         'purity_completeness_per_bin']


def load_summary(input_dir):
    f = open(input_dir + "/summary.tsv", 'r')
    field_names = f.readline().rstrip('\n').split('\t')
    table = []
    for line in f:
        values = line.rstrip('\n').split('\t')
        table.append(dict(zip(field_names, values)))
    f.close()
    return table


def create_completeness_contamination_table(input_dir, output_dir, gold_standard=None):
    summary = load_summary(input_dir)
    with open(output_dir + "/completeness_contamination_table.tex", "w") as texfile:
        texfile.write('{}{}'.format(r'\documentclass{article}', '\n'))
        texfile.write('{}{}'.format(r'\usepackage[T1]{fontenc}', '\n'))
        texfile.write('{}{}'.format(r'\usepackage{multirow}', '\n'))
        texfile.write('{}{}'.format(r'\pagestyle{empty}', '\n'))
        texfile.write('{}{}'.format(r'\begin{document}', '\n'))
        texfile.write('{}{}'.format(r'\fontfamily{cmss}\selectfont', '\n'))
        texfile.write('{}{}'.format(r'\begin{tabular}{l|r|r|r|r}', '\n'))
        texfile.write('{}{}'.format(r'\multicolumn{2}{c|}{\textbf{Genome binner}}&\multicolumn{3}{|c}{\textbf{Predicted bins}}\\', '\n'))
        texfile.write('{}{}'.format(r'\multicolumn{2}{c|}{\textbf{(\% contamination)}}&\multicolumn{3}{|c}{\textbf{(\% completeness)}}\\', '\n'))
        texfile.write('{}{}'.format(r'\cline{3-5}', '\n'))
        texfile.write('{}{}'.format(r'\multicolumn{2}{c|}{}&>50\%&>70\%&>90\%\\', '\n'))
        texfile.write('{}{}'.format(r'\hline', '\n'))

        if gold_standard:
            num_genomes = str(len(gold_standard.genome_id_to_list_of_contigs))
            texfile.write('{}{}'.format(r'\multicolumn{2}{l|}{}&&&\\[-1em]', '\n'))
            texfile.write('{}{}'.format(r'\multicolumn{2}{l|}{Gold standard}&' + num_genomes + r'&' + num_genomes + r'&' + num_genomes + r'\\', '\n'))
            texfile.write('{}{}'.format(r'\hline', '\n'))

        for row in summary:
            texfile.write('{}{}'.format(r'\multirow{2}{*}{' + row['tool'] + r'}&<10\%&' + row['>0.5compl<0.1cont'] + r'&' + row['>0.7compl<0.1cont'] + r'&' + row['>0.9compl<0.1cont'] + r'\\', '\n'))
            texfile.write('{}{}'.format(r'&<5\%&' + row['>0.5compl<0.05cont'] + r'&' + row['>0.7compl<0.05cont'] + r'&' + row['>0.9compl<0.05cont'] + r'\\', '\n'))
            texfile.write('{}{}'.format(r'\hline', '\n'))
        texfile.write('{}{}'.format(r'\end{tabular}', '\n'))
        texfile.write('{}{}'.format(r'\end{document}', '\n'))
    texfile.close()
    os.system("latex -output-directory " + output_dir + " " + output_dir + "/completeness_contamination_table.tex")
    os.system("dvips " + output_dir + "/completeness_contamination_table.dvi " + "-o " + output_dir + "/completeness_contamination_table.ps")
    os.system("ps2eps -f " + output_dir + "/completeness_contamination_table.ps")
    os.system("cp " + output_dir + "/completeness_contamination_table.eps " + input_dir)


def create_pdf(input_dir, output_dir, plots):
    with open(output_dir + "/summary.tex", "w") as texfile:
        texfile.write('{}{}'.format(r'\documentclass{article}', '\n'))
        texfile.write('{}{}'.format(r'\usepackage[left=0.3in,right=0.3in,top=0.5in,bottom=0.3in]{geometry}', '\n'))
        texfile.write('{}{}'.format(r'\usepackage{graphicx}', '\n'))
        texfile.write('{}{}'.format(r'\usepackage{floatrow}', '\n'))
        texfile.write('{}{}'.format(r'\usepackage{subfig}', '\n'))
        texfile.write('{}{}'.format(r'\usepackage[T1]{fontenc}', '\n'))
        texfile.write('{}{}'.format(r'\usepackage{multirow}', '\n'))
        texfile.write('{}{}'.format(r'\floatsetup[figure]{style=plain,subcapbesideposition=top}', '\n'))
        texfile.write('{}{}'.format(r'\begin{document}', '\n'))
        texfile.write('{}{}'.format(r'\begin{figure}[!ht]', '\n'))
        texfile.write('{}{}'.format(r'\centering', '\n'))
        texfile.write('{}{}'.format(r'\includegraphics[scale=0.75]{' + input_dir + '/legend}\\', '\n'))
        texfile.write('{}{}'.format(r'\ffigbox', '\n'))
        texfile.write('{}{}'.format(r'{', '\n'))
        i = 1
        for plot in plots:
            if i % 2 != 0:
                texfile.write('{}{}'.format(r'\begin{subfloatrow}', '\n'))
            texfile.write('{}{}'.format(r'\sidesubfloat[]{\includegraphics[width=0.45\linewidth]{' + input_dir + '/' + PLOTS[plot - 1] + r'}}', '\n'))
            if i % 2 == 0 or i == len(plots):
                texfile.write('{}{}'.format(r'\end{subfloatrow}', '\n'))
                if i < len(plots):
                    texfile.write('{}{}'.format(r'\par\bigskip', '\n'))
            i += 1
        texfile.write('{}{}'.format(r'}', '\n'))
        texfile.write('{}{}'.format(r'{}', '\n'))
        texfile.write('{}{}'.format(r'\end{figure}', '\n'))
        texfile.write('{}{}'.format(r'\end{document}', '\n'))
    texfile.close()
    os.system("latex -output-directory " + output_dir + " " + output_dir + "/summary.tex")
    os.system("dvipdf " + output_dir  + "/summary.dvi " + output_dir + "/summary.pdf")


def main():
    if os.name == 'nt':
        exit('Sorry, this tool does not support Windows.')
    parser = argparse.ArgumentParser(description="Combine figures and table of completeness and contamination into file summary.pdf")
    parser.add_argument('-c', '--figure_codes', help="Comma-separated figure codes without spaces (1=average purity vs. average completeness, 2=ari vs. %%assigned bps, 3=weighed purity vs. weighed completeness, 4=table of contamination and completeness, 5=purity vs. completeness per bin; default 1,2,3,4)", required=False)
    parser.add_argument('-i', '--input_dir', help="Directory to read file summary.tsv and figures from", required=True)
    parser.add_argument('-o', '--output_dir', help="Directory to write pdf and temporary files (if different from input_dir)", required=False)
    parser.add_argument('-g', '--gold_standard_file', help=argparse_parents.HELP_GOLD_STANDARD_FILE, required=False)
    parser.add_argument('-r', '--remove_genomes', help=argparse_parents.HELP_GENOMES_FILE, required=False)
    parser.add_argument('-k', '--keyword', help=argparse_parents.HELP_KEYWORD, required=False)

    args = parser.parse_args()
    figure_codes = [1, 2, 3, 4]
    output_dir = os.path.abspath(args.input_dir)
    input_dir = os.path.abspath(args.input_dir)

    if args.figure_codes:
        figure_codes = [int(x.strip()) for x in args.figure_codes.split(',')]
    if args.output_dir:
        output_dir = os.path.abspath(args.output_dir)
        load_data.make_sure_path_exists(args.output_dir)

    gold_standard = None
    if args.gold_standard_file:
        gold_standard = load_data.get_genome_mapping_without_lenghts(args.gold_standard_file, args.remove_genomes, args.keyword)

    if 4 in figure_codes:
        create_completeness_contamination_table(input_dir, output_dir, gold_standard)

    create_pdf(input_dir, output_dir, figure_codes)


if __name__ == "__main__":
    main()
