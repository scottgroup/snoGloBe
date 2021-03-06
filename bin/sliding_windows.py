#!/usr/bin/env python3

from Bio import SeqIO
import pandas as pd
import sys


def read_trx_fasta(fasta):
    trx_dict = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
    # set start and strand decoy values for compatibility with target_dict
    trx_dict = {k: {
        'seq': str(v.seq).upper().replace('T', 'U'),
        'start' : 0,
        'strand': '+'
    } for k, v in trx_dict.items() if len(v.seq) > 13}
    return trx_dict


def seq_onehotmatrix(df):
    df_seq = df.seq.str.extractall('(.)')[0].unstack()
    df_seq = pd.get_dummies(df_seq)
    df = df.drop(['seq'], axis=1)
    df = df.sort_index()
    df_seq = df_seq.sort_index()
    df = df.join(df_seq)
    return df


def sliding_window(seq_dict, step):
    length = 13
    windows_dict = {}
    for trx, v in seq_dict.items():
        for i in range(0, len(v['seq']) - length + 1, step):
            if v['strand'] == '+':
                window_start = v['start'] + i
            elif v['strand'] == '-':
                window_start = v['start'] + len(v['seq']) - length - i
            else:
                print('Error! Wrong choice of strand : %s\n'
                      'Acceptable choices are : + or -' % v['strand'], file=sys.stderr)
                sys.exit(1)
            windows_dict[trx + '_' + str(window_start)] = [
                str(v['seq'][i:i + length]).upper(),
                i,
                len(v['seq'])
            ]
    df = pd.DataFrame.from_dict(windows_dict, orient='index', columns=['seq', 'idx', 'length'])
    df['seq'] = df.seq.str.upper()
    df['seq'] = df.seq.str.replace('T', 'U')
    df = seq_onehotmatrix(df)
    return df
