#!/usr/bin/env python3

import os
import sys
from distutils.spawn import find_executable


def check_gtf(gtf_file):
    if gtf_file.rsplit('.',1)[-1] != 'gtf':
        print('Error: bad annotation file extension. '
              'Please use a gene transfer format (.gtf) as input annotation.',
              file=sys.stderr)
        sys.exit(1)
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#') == False:
                if line.split('\t')[2] == 'gene':
                    gene_test_line = line
                elif line.split('\t')[2] == 'transcript':
                    transcript_test_line = line
                elif line.split('\t')[2] == 'exon':
                    exon_test_line = line
                try:
                    gene_test_line
                    transcript_test_line
                    exon_test_line
                except:
                    pass
                else:
                    break
    # Check if any gene, transcript and exon entry was found in the gtf file.
    try:
        gene_test_line
    except:
        print(
            'Error: bad annotation format, no gene entry.'
            ' Please use a gene transfer format (.gtf) as input annotation.',
            file=sys.stderr)
        sys.exit(1)
    try:
        transcript_test_line
    except:
        print(
            'Error: bad annotation format, no transcript entry. '
            'Please use a gene transfer format (.gtf) as input annotation.',
            file=sys.stderr)
        sys.exit(1)
    try:
        exon_test_line
    except:
        print(
            'Error: bad annotation format, no exon entry. '
            'Please use a gene transfer format (.gtf) as input annotation.',
            file=sys.stderr)
        sys.exit(1)

    if 'gene_biotype' not in gene_test_line:
        print('Error: There is no "gene_biotype" entry for genes.'
              ' Please use a gene transfer format (.gtf) annotation file obtained from Ensembl to correct this.',
              file=sys.stderr)
        sys.exit(1)
    if 'gene_id' not in gene_test_line:
        print('Error: There is no "gene_id" entry for genes.'
              ' Please use a gene transfer format (.gtf) annotation file obtained from Ensembl to correct this.',
              file=sys.stderr)
        sys.exit(1)


def check_output(output):
    if os.path.isdir(output):
        print('Error: output specified is a directory, please provide a file name to be created',
              file=sys.stderr)
        sys.exit(1)
    output_path = os.path.dirname(output)
    cwd = os.getcwd()
    if os.path.isdir(output_path) == False:
        try:
            if os.path.isdir(output_path) == False:
                if os.path.isdir(os.path.join(cwd, output_path)) == False:
                    sys.exit(1)
        except:
            print("Error: path to output does not exist: %s. Please specify a valid output path." % (output),
                  file=sys.stderr)
            sys.exit(1)


def check_dependencies():
    exist = find_executable('bedtools')
    if exist == None:
        print('Error: bedtools is not installed. Please read the README.md file '
              'for more information about snoGloBe\'s prerequisites.',
              file=sys.stderr)
        sys.exit(1)
