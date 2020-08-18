#!/usr/bin/env python3

__author__ = "Gabrielle Deschamps-Francoeur"
__email__ = "gabrielle.deschamps-francoeur@usherbrooke.ca"
__version__ = '0.1'

import argparse
from sliding_windows import sliding_window, read_trx_fasta
from make_feature_bed import make_bed
import os
import pandas as pd
import extract_features as ef
from predict_mp import make_pred
from consecutive_windows import cons_windows
import sys
import tests


def read_gtf(gtf_file):
    df_gtf = pd.read_csv(gtf_file, sep='\t', comment='#',
                         names=['seqname', 'source', 'feature', 'start', 'end', 'score',
                                'strand', 'frame', 'attributes'],
                         dtype={'seqname': str, 'start': int, 'end': int},
                         usecols=['seqname', 'feature', 'start', 'end', 'score', 'strand', 'attributes'])
    df_gtf['gene_id'] = df_gtf['attributes'].str.extract('gene_id\s"([^;]+);?"', expand=True)
    df_gtf['gene_biotype'] = df_gtf['attributes'].str.extract('gene_biotype\s"([^;]+);?"', expand=True)
    df_gtf['transcript_id'] = df_gtf['attributes'].str.extract('transcript_id\s"([^;]+);?"', expand=True)
    df_gtf = df_gtf.drop(['attributes'], axis=1)
    return df_gtf


def sno_windows(sno_fasta, sno_file):
    # read sno fasta file and split in sliding windows
    df_sno = sliding_window(sno_fasta, 1)
    df_sno['rel_pos'] = df_sno.idx / df_sno.length
    df_sno['rel_pos'] = df_sno.rel_pos.round(3)
    df_sno = df_sno.drop(['idx', 'length'], axis=1)
    df_sno['snoid'] = df_sno.index
    if len(df_sno[df_sno.snoid.str.contains(',')]) > 0:
        print('Warning: commas (,) in snoRNA ids will be changed for dot (.)')
        df_sno['snoid'] = df_sno['snoid'].str.replace(',','.')
    df_sno = df_sno[df_sno.columns.sort_values()]
    df_sno.to_csv(sno_file, index=False)


def target_windows(target_fasta, df_gtf, step, target_file, bedfile):
    df_target = sliding_window(target_fasta, step)
    df_target['gene_id'] = df_target.index.str.rsplit('_', 1).str[0]
    df_gene = df_gtf[df_gtf.feature == 'gene'][['gene_id', 'seqname', 'start', 'end', 'strand']]
    df_target = df_target.reset_index().merge(df_gene, on='gene_id', how='left').set_index('index')
    df_target.loc[df_target.strand == '+', 'window_start'] = df_target.start + df_target.idx
    df_target.loc[df_target.strand == '+', 'window_end'] = df_target.window_start + 13
    df_target.loc[df_target.strand == '-', 'window_start'] = df_target.end - df_target.idx - 13 + 1
    df_target.loc[df_target.strand == '-', 'window_end'] = df_target.end - df_target.idx + 1
    df_target = df_target.drop(['start', 'end', 'gene_id', 'idx', 'length'], axis=1)
    df_bed = df_target[['seqname', 'window_start', 'window_end', 'strand']].copy(deep=True)
    df_target = df_target.drop(['seqname', 'strand', 'window_start', 'window_end'], axis=1)
    df_bed['score'] = 0
    df_bed['idx'] = df_bed.index
    df_out = ef.extract_features(df_bed, os.path.dirname(target_file), bedfile)
    df_out = df_target.merge(df_out, left_index=True, right_index=True)
    cols = ef.column_list(df_gtf)
    df_out = ef.set_columns(df_out, cols)
    df_out['wid'] = df_out.index
    df_out.to_csv(target_file, index=False)


def extract_seq(row, dict, seq_type):
    # Sequences in 5'-3' orientation
    seq = dict[row[seq_type + '_id']].seq
    start = row[seq_type + '_window_start']
    end = row[seq_type + '_window_end']
    window_seq = seq[start:end]
    return window_seq


def interaction_sequence(outfile, sno_fasta, target_fasta, chunksize):
    tempfile = outfile + '.seq'
    target_dict = read_trx_fasta(target_fasta)
    sno_dict = read_trx_fasta(sno_fasta)
    wlength = 13
    for i, df in enumerate(pd.read_csv(outfile, chunksize=chunksize)):
        if 'target_window_end' not in df.columns:
            df['target_window_end'] = df['target_window_start'] + wlength
            df['sno_window_end'] = df['sno_window_start'] + wlength
        if i == 0:
            mode = 'w'
            header = True
        else:
            mode = 'a'
            header = False
        if df.empty is False:
            df['target_seq'] = df.apply(lambda row: extract_seq(row, target_dict, 'target'), axis=1)
            df['sno_seq'] = df.apply(lambda row: extract_seq(row, sno_dict, 'sno'), axis=1)
            score_cols = [c for c in df.columns if 'score' in c]
            df[score_cols] = df[score_cols].round(3)
            if len(score_cols) > 1:
                score_cols = ['count'] + score_cols
            sorted_cols = ['target_id', 'target_window_start', 'target_window_end', 'sno_id', 'sno_window_start',
                           'sno_window_end'] + score_cols + ['target_seq', 'sno_seq']
            df = df[sorted_cols]
            df.to_csv(tempfile, index=False, mode=mode, header=header)
    os.rename(tempfile, outfile)


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("sno_fasta", help="fasta file containing snoRNA sequences")
    parser.add_argument("target_fasta", help="fasta file containing target gene sequences, "
                                             "gene_id should match the ones in the gtf file")
    parser.add_argument("gtf", help="Annotation file in .gtf format. Preferably an annotation of the whole "
                                    "genome or whole chromosomes of specified targets")
    parser.add_argument("output", help="Output file name")
    parser.add_argument("-v", "--version", help="Show the version and exit",
                        action="store_true")
    parser.add_argument("-n", "--nb_threads", help="Number of threads to be used. Default: 1",
                        type=int, default=1)
    parser.add_argument("-s", "--stepsize", help="Number of nucleotides between each target sliding window. Default: 1",
                        type=int, default=1)
    parser.add_argument("-c", "--chunksize", help="Number of combinations to run at once. Default: 3000000",
                        type=int, default=3000000)
    parser.add_argument("-t", "--threshold", help="Minimum score for an interaction to be reported, "
                                                  "value between 0 and 1. Default: 0.95",
                        type=float, default=0.95)
    parser.add_argument("-m", "--merge", help="Use this option to merge consecutive positive windows as one",
                        action='store_true')
    parser.add_argument("-w", "--nb_windows", help="Minimum number of consecutive windows for an "
                                                   "interaction to be reported. Must be used with -m option. Default: 1",
                        type=int, default=1)
    parser.add_argument("--seq", help="Add target and snoRNA interaction sequences to the output",
                        action='store_true')

    if sys.argv[1] in ['-v', '--version']:
        print('snoGloBe version :', __version__)
        exit(0)

    args = parser.parse_args()

    sno_fasta = args.sno_fasta
    target_fasta = args.target_fasta
    gtf = args.gtf
    nb_threads = args.nb_threads
    step = args.stepsize
    chunksize = args.chunksize
    threshold = args.threshold
    outfile = os.path.abspath(args.output)
    outpath = os.path.dirname(outfile)
    merge_conswindows = args.merge
    nb_windows = args.nb_windows
    add_seq = args.seq

    tests.check_dependencies()
    tests.check_gtf(gtf)
    tests.check_output(outfile)

    if nb_windows != 1 and merge_conswindows is False:
        print('-w/--nb_windows must be used with -m/--merge option. -w/--windows %d will be ignored' % nb_windows)

    if threshold < 0 or threshold >1:
        print("Wrong choice of -t/--threshold. Value must be between 0 and 1, chosen value was: ", threshold,
              file=sys.stderr)
        exit(1)

    # define temp file names
    sno_file = os.path.join(outpath, os.path.basename(sno_fasta).replace('.fa', '.input.csv'))
    bedfile = os.path.join(outpath, os.path.basename(gtf).replace('.gtf', '.features.bed'))
    target_file = os.path.join(outpath, os.path.basename(target_fasta).replace('.fa', '.input.csv'))

    df_gtf = read_gtf(gtf)
    sno_windows(sno_fasta, sno_file)
    make_bed(df_gtf, bedfile)
    target_windows(target_fasta, df_gtf, step, target_file, bedfile)
    make_pred(step, sno_file, target_file, chunksize, nb_threads, outfile, threshold, merge_conswindows)

    # remove temp files
    os.remove(sno_file)
    os.remove(target_file)
    os.remove(bedfile)

    # merge consecutive windows
    if merge_conswindows is True:
        cons_windows(outfile, nb_windows, chunksize, step, nb_threads)

    if add_seq is True:
        interaction_sequence(outfile, sno_fasta, target_fasta, chunksize)


if __name__ == '__main__':
    main()
