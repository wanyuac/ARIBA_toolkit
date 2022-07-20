#!/usr/bin/env python

"""
This script extracts assembled gene sequences from input FASTA files and pools them into single FASTA files
according to gene names. Note that, a 'gene' corresponds to a sequence 'cluster' in ARIBA's reference database.

Input filename: [isolate name]__genes.fna ('__' is essential, whereas 'genes.fna' is replaceable).

Example command:
    python parseGeneSeqs.py --input ariba_output/*__genes.fna --outdir alleles > sequence_summary.tsv

Copyright (C) 2022 Yu Wan <wanyuac@126.com>
Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
Publication: 20 July 2022; last modification: 20 July 2022
"""

import os
from Bio import SeqIO
from argparse import ArgumentParser

def parse_arguments():
    parser = ArgumentParser(description = "Extract assembled gene sequences from ARIBA's outputs")
    parser.add_argument('-i', '--input', dest = 'input', nargs = '+', type = str, required = True, help = "Input FASTA files")
    parser.add_argument('-o', '--outdir', dest = 'outdir', type = str, default = 'output', required = False, help = "Output directory")
    
    return parser.parse_args()


def main():
    args = parse_arguments()
    
    # Output settings
    outdir = args.outdir
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    print('\t'.join(['Isolate', 'Gene', 'Length']))  # The header line of the output summary table
    
    # Processes each input FASTA file
    records = dict()
    for f in args.input:
        f_basename = os.path.basename(f)
        i, _ = f_basename.split('__')  # Get the isolate name
        for seq in SeqIO.parse(f, 'fasta'):
            gene = seq.id.split('.')[0]
            seq.id = gene + '__' + i  # Append the isolate name to the sequence name
            print('\t'.join([i, gene, str(len(str(seq.seq)))]))
            if gene in records.keys():
                records[gene].append(seq)
            else:  # A new gene is encountered
                records[gene] = [seq]
    
    # Writes sequences
    for g, seqs in records.items():
        SeqIO.write(seqs, os.path.join(outdir, g + '.fna'), 'fasta')

    return

if __name__ == '__main__':
    main()
