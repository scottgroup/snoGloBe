#!/usr/bin/env python3

from Bio import SeqIO
import pandas as pd


def read_trx_fasta(fasta):
    trx_dict = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
    trx_dict = {k: v for k, v in trx_dict.items() if len(v.seq) > 13}
    return trx_dict


def seq_onehotmatrix(df):
    df_seq = df.seq.str.extractall('(.)')[0].unstack()
    df_seq = pd.get_dummies(df_seq)
    df = df.drop(['seq'], axis=1)
    df = df.merge(df_seq, right_index=True, left_index=True)
    return df


def sliding_window(fasta, step):
    seq_dict = read_trx_fasta(fasta)
    length = 13
    windows_dict = {}
    for trx, v in seq_dict.items():
        for i in range(0, len(v.seq) - length + 1, step):
            windows_dict[trx + '_' + str(i)] = [
                str(v.seq[i:i + length]).upper(),
                i,
                len(v.seq)
            ]
    df = pd.DataFrame.from_dict(windows_dict, orient='index', columns=['seq', 'idx', 'length'])
    df['seq'] = df.seq.str.upper()
    df['seq'] = df.seq.str.replace('T', 'U')
    df = seq_onehotmatrix(df)
    return df
