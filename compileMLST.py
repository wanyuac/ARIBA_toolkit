"""
Compile results of MLST analysis using ARIBA.

This script assumes subdirectory are named after genomes. It takes as input from stdin and
produces output files under the current working directory.

Command line:
    find ./* -maxdepth 1 -type d | python compileMLST.py
    cat outdir.txt | python compileMLST.py -p mlst
    
Note that do not use `find . -maxdepth 1 -type d` or `find ./ -maxdepth 1 -type d` for this
script because both commands return the current working directory './' in the first line.

Python 3 is preferred to run this script, though it is compatible to Python 2.

Copyright 2020 Yu Wan <wanyuac@126.com>
Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
Publication: 8 Mar 2020
"""

from __future__ import print_function
import os
import sys
from argparse import ArgumentParser


def parse_arguments():
    parser = ArgumentParser(description= 'Compile ARIBA outputs of MLST analysis')
    parser.add_argument('-p', '--prefix', dest = 'p', type = str, required = False, default = '',\
                        help = 'Prefix of output files')
    
    return parser.parse_args()


def main():
    # Define types of files to be processed
    FILES = {'STs' : 'mlst_report.tsv', \
             'hits' : 'mlst_report.details.tsv', \
             'scores' : 'report.tsv'}
    
    args = parse_arguments()
    genomes = get_genome_names(sys.stdin)
    
    for file_type, file_name in FILES.items():
        compile_files(file_type, file_name, genomes, args.p)
        
    return


def get_genome_names(paths):
    """
    Extracts genome names from paths of directories containing ARIBA's output
    """
    genomes = {}
    
    for p in paths:
        p = p.rstrip('\n')  # Otherwise, every genome name has a newline character.
        g = os.path.basename(p)  # For example, './genome1' becomes 'genome1'.
        if g == '':
            print(p + ' is not parsed as it points to the current working directory.')
        else:
            genomes[g] = p
            
    if len(genomes) == 0:
        sys.exit('Error: No genome name is extracted from paths.')
    
    return genomes


def compile_files(report_type, report_name, genomes, prefix):
    """
    Compiles files of a given type.
    """
    print('Start to compile report type ' + report_type + '.')
    file_name = report_type + '.tsv'
    if prefix != '':
        file_name = prefix + '__' + file_name
    compiled_report = open(file_name, 'w')  # Output file
    
    report_count = 0  # File counter
    for g, p in genomes.items():  # g: genome ID; p: path to a directory
        try:
            report = open(os.path.join(p, report_name), 'r')
            report_count += 1
            lines = report.readlines()  # Does not need to use read().splitlines() to strip off newline characters.
            line_count = len(lines)
            
            if report_count == 1:  # Read the header of the first file, assuming other files have the same format
                header = lines[0]
                if header.startswith('#'):
                    header = header[1 : len(header)]  # Chop away the leading hash character.
                compiled_report.write('genome' + '\t' + header)
            
            if line_count > 1:  # Skip the header line for the block to be written into the output file
                lines = lines[1 : line_count]
                for r in lines:  # Write into the output file
                    compiled_report.write(g + '\t' + r)
            else:
                print('Warning: file ' + p + 'only has one line. Hence this file is skipped.')
            
            report.close()
        except:
            compiled_report.close()
            sys.exit('Error: Input file ' + p + ' or output file ' + file_name + ' is not accessible.')
    
    compiled_report.close()
    
    return


if __name__ == '__main__':
    main()
