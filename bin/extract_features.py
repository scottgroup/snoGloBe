#!/usr/bin/env python3

import os
import subprocess
import sys


def column_list():
    cols = ['Mt_rRNA', 'Mt_tRNA', 'TEC', 'TR_V_gene', 'exon', 'five_prime_utr',
            'intron', 'lncRNA', 'miRNA', 'misc_RNA', 'protein_coding', 'pseudogene',
            'rRNA', 'ribozyme', 'sRNA', 'scRNA', 'scaRNA', 'snRNA', 'snoRNA', 'tRNA',
            'three_prime_utr', 'vaultRNA', 'wid']
    length = 13
    nts = ['A', 'U', 'C', 'G']
    seq_cols = [str(i) + '_' + str(j) for i in range(0, int(length)) for j in nts]
    cols.extend(seq_cols)
    cols.sort()
    return cols


def set_columns(df):
    cols = column_list()
    missing_cols = set(cols) - set(df.columns)
    unexpected_cols = set(df.columns) - set(cols)
    if len(unexpected_cols) > 0:
        print('Warning: The following features won\'t be taking into account for the prediction. '
              'Please refer to the manual for a list of accepted features.', file=sys.stderr)
        print(unexpected_cols, file=sys.stderr)
    for i in missing_cols:
        df[i] = 0
    df = df[cols]
    new_cols = {i: str(int(i.split('_')[0]) + 13) + '_' + '_'.join(i.split('_')[1:])
    if i.split('_')[0].isdigit() else i for i in cols}
    df = df.rename(columns=new_cols)
    return df


def extract_features(fpath, bedfile, target_file):
    temp_bed = fpath + '.targetwindows.seq.bed'

    bedtools_cmd = ['bedtools', 'intersect', '-a', temp_bed,
                    '-b', bedfile, '-s', '-wb']
    cols = ['Mt_rRNA', 'Mt_tRNA', 'TEC', 'TR_V_gene', 'exon', 'five_prime_utr',
            'intron', 'lncRNA', 'miRNA', 'misc_RNA', 'protein_coding', 'pseudogene',
            'rRNA', 'ribozyme', 'sRNA', 'scRNA', 'scaRNA', 'snRNA', 'snoRNA', 'tRNA',
            'three_prime_utr', 'vaultRNA']

    cut_cmd = ('cut', '-f', '4,7,14-' + str(14 + len(cols) - 1))

    groupby_cmd = ('bedtools', 'groupby',
                   '-g', '1,2',
                   '-c', ','.join(str(i) for i in range(3, 3 + len(cols))),
                   '-o', 'max')

    res0 = subprocess.Popen(bedtools_cmd, stdout=subprocess.PIPE)
    res1 = subprocess.Popen(cut_cmd, stdin=res0.stdout, stdout=subprocess.PIPE)
    res2 = subprocess.Popen(groupby_cmd, stdin=res1.stdout, stdout=open(target_file, mode='a'))
    res2.wait()
    res1.wait()
    res0.wait()
    if res0.returncode == 0 and res1.returncode == 0 and  res2.returncode == 0:
        os.remove(temp_bed)
    else:
        sys.exit(1)
