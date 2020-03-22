#!/usr/bin/env python

"""
Concatenate allele sequences of a given list of STs. It is useful for comparing STs at the nucleotide
level and MLST-based phylogenetic construction. Currently this script does not support allele variants,
namely, requires exact hits of query genomes to reference allele sequences.

Example command lines:
    python concMLSTalleles.py -p st_profiles.tsv -s ./sequence_dir -e tfa > mlstSeqs.fna
    python concMLSTalleles.py -p st_profiles_new.tsv -s ./sequence_dir -e tfa >> mlstSeqs.fna  # Append sequences

This script assumes input ST profiles are offered in a tab-delimited file with columns of STs and
names of MLST loci. (No other columns should be allowed) For example:

ST    locus_1    locus_2    locus_3    locus_4    locus_5    locus_6    locus_7
1    11    21    31    41    51    61    71
2    12    21    31    41    51    61    71
...

The directory of MLST allele sequences is expected to be the same as that downloaded from the PubMLST
database using ARIBA. (pubmlst_download). The '-s' parameter does not require an end forward slash
for the path. By default, the delimiter for locus name (for instance, gene adk) and allele number in
an allele name (for example, adk.1) is a period, although some databases use an underscore as the
delimiter.

Dependencies: Python v3, packages pandas (pandas.pydata.org/) and biopython.

Copyright 2020 Yu Wan <wanyuac@126.com>
Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
Publication: 22 Mar 2020
"""

import os
import sys
import pandas as pd
from argparse import ArgumentParser
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import single_letter_alphabet
from collections import defaultdict


def parse_arguments():
    parser = ArgumentParser(description= 'Concatenate MLST alleles according to ST profiles')
    parser.add_argument('-p', dest = 'p', type = str, required = True,\
                        help = 'ST profiles of interest')
    parser.add_argument('-s', dest = 's', type = str, required = True,\
                        help = 'Path to the directory containing allele sequences of every MLST locus')
    parser.add_argument('-d', dest = 'd', type = str, required = False, default = '.',\
                        help = 'Common delimiter between the locus name and allele number in an allele name')
    parser.add_argument('-e', dest = 'e', type = str, required = False, default = 'tfa',\
                        help = 'Filename extension of input FASTA files. (tfa, fasta, fna, etc)')
    parser.add_argument('-v', dest = 'v', action = 'store_true', required = False,\
                        help = 'Flag this option to print sequence headers in a verbose way (with allele information)')

    return parser.parse_args()


def main():
    args = parse_arguments()
    
    # Import ST profiles of interest as a data frame (type: pandas.core.frame.DataFrame)
    try:
        st_prof = pd.read_csv(args.p, sep = '\t', header = 0, index_col = 0, dtype = str)  # Note that the first row and column start from index zero.
        st_prof.index = st_prof.index.map(str)  # Map pandas.core.indexes.numeric.Int64Index to strings.
    except:
        sys.exit("Error: cannot read ST profiles.")
    
    # Import MLST allele sequences
    loci = list(st_prof.columns)
    db = import_mlst_alleles(args.s, loci, args.e)  # Read FASTA files of MLST alleles
    
    # Concatenate and print allele sequences of every ST of interest
    sts = list(st_prof.index)  # Row names: STs
    for st in sts:
        seq_record = conc_seqs(st_prof.loc[st, : ], loci, db, args.d, args.v)  # Feed a row of ST profiles into the function
        if args.v:
            print('>%s %s %s\n%s' % (seq_record.id, seq_record.name, seq_record.description, seq_record.seq))  # Print this SeqRecord
        else:
            print('>%s\n%s' % (seq_record.id, seq_record.seq))  # Print sequence headers in a concise manner (Compatible to SplitTree4)
    
    return


def import_mlst_alleles(seq_dir, loci, ext):
    """
    Import FASTA files of MLST allele sequences.
    Output data structure: {locus : {allele : sequence}}
    """
    db = defaultdict(dict)
    for locus in loci:
        fasta = os.path.join(seq_dir, '.'.join([locus, ext]))
        if os.path.exists(fasta):
            db[locus] = SeqIO.to_dict(SeqIO.parse(fasta, 'fasta'))  # Every element dictionary is a SeqRecord object.
        else:
            sys.exit("Error: sequence file %s is not accessible." % fasta)
    
    return db


def conc_seqs(pf, loci, db, delim, verbose):
    """
    Concatenate allele sequences for a given ST.
    """
    seqs = []  # A list of Seq objects
    descr = []
    for locus in loci:
        allele = locus + delim + pf[locus]  # For example, 'locus_1'.
        record = db[locus][allele]  # Record is a SeqRecord object.
        seqs.append(record.seq)
        descr.append(allele + '(%i)' % len(record.seq))  # locus_1(length)
    seq = Seq('', single_letter_alphabet).join(seqs)  # Concatenate allele sequences
    if verbose:
        seq_record = SeqRecord(seq, id = pf.name, name = str(len(seq)) + 'bp', description = ','.join(descr))
    else:
        seq_record = SeqRecord(seq, id = pf.name, name = '', description = '')
    
    return seq_record


if __name__ == '__main__':
    main()
