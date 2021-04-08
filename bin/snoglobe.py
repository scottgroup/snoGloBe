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
from fetch_sequence import get_sequence_from_coordinates as get_seq
import re


def read_gtf(gtf_file, verbose):
    if verbose:
        print('Reading gtf')
    gtf_list = []
    with open(gtf_file) as f:
        for line in f:
            if line[0] == '#':
                continue
            line = line.split('\t')
            if line[2] not in ['gene', 'transcript', 'exon', 'three_prime_utr', 'five_prime_utr']:
                continue
            cols_to_keep = [0, 2, 3, 4, 5, 6]
            attribute = line[8]
            gtf_entry = [item for idx, item in enumerate(line) if idx in cols_to_keep]
            gene_id = re.search('gene_id\s"([^;]+);?"', attribute).group(1)
            try:
                transcript_id = re.search('transcript_id\s"([^;]+);?"', attribute).group(1)
            except AttributeError:
                transcript_id = None
            try:
                exon_id = re.search('exon_id\s"([^;]+);?"', attribute).group(1)
            except AttributeError:
                exon_id = None
            try:
                gene_biotype = re.search('gene_biotype\s"([^;]+);?"', attribute).group(1)
            except AttributeError:
                gene_biotype = None
            gtf_entry.extend([gene_id, transcript_id, exon_id, gene_biotype])
            gtf_list.append(gtf_entry)
    df_gtf = pd.DataFrame(gtf_list,
                          columns=['seqname', 'feature', 'start', 'end', 'score', 'strand',
                                   'gene_id', 'transcript_id', 'exon_id', 'gene_biotype'])
    df_gtf[['start', 'end']] = df_gtf[['start', 'end']].astype(int)
    df_gtf.loc[df_gtf.feature == 'gene', 'feature_id'] = df_gtf.gene_id
    df_gtf.loc[df_gtf.feature == 'transcript', 'feature_id'] = df_gtf.transcript_id
    df_gtf.loc[df_gtf.feature == 'exon', 'feature_id'] = df_gtf.exon_id
    return df_gtf


def sno_windows(sno_fasta, sno_file, verbose):
    # read sno fasta file and split in sliding windows
    if verbose:
        print('preparing snoRNA windows')
    sno_dict = read_trx_fasta(sno_fasta)
    df_sno = sliding_window(sno_dict, 1)
    df_sno['rel_pos'] = df_sno.idx / df_sno.length
    df_sno['rel_pos'] = df_sno.rel_pos.round(3)
    df_sno = df_sno.drop(['idx', 'length'], axis=1)
    df_sno['snoid'] = df_sno.index
    if len(df_sno[df_sno.snoid.str.contains(',')]) > 0:
        print('Warning: commas (,) in snoRNA ids will be changed for dot (.)', file=sys.stderr)
        df_sno['snoid'] = df_sno['snoid'].str.replace(',','.')
    df_sno = df_sno[df_sno.columns.sort_values()]
    df_sno.to_csv(sno_file, index=False)


def target_windows(target_dict, df_gtf, step, target_file, bedfile):
    df_target = sliding_window(target_dict, step)
    df_target['feature_id'] = df_target.index
    df_target[['feature_id', 'seqname', 'strand', 'window_start']] = df_target.feature_id.str.rsplit('_', 3, expand=True)
    df_target['window_end'] = df_target.window_start.astype(int) + 13
    df_target = df_target.drop(['idx', 'length', 'feature_id'], axis=1)
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
    feature_id = row[seq_type + '_id']
    strand = dict[feature_id]['strand']
    seq = dict[feature_id]['seq']
    gene_start = dict[feature_id]['start']
    gene_end = gene_start + len(seq)
    if strand == '+' :
        gene_start = dict[feature_id]['start']
        start = row[seq_type + '_window_start'] - gene_start # relative start
        end = row[seq_type + '_window_end'] - gene_start  # relative end
    elif strand == '-':
        start = gene_end - row[seq_type + '_window_end'] # relative start
        end = gene_end - row[seq_type + '_window_start'] # relative end
    else:
        print('wrong choice of strand: %s' % strand, file=sys.stderr)
        sys.exit(1)
    window_seq = seq[start:end]
    return window_seq


def interaction_sequence(outfile, sno_fasta, target_dict, chunksize):
    tempfile = outfile + '.seq'
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
        else:
            df['target_seq'] = None
            df['sno_seq'] = None
        score_cols = [c for c in df.columns if 'score' in c]
        df[score_cols] = df[score_cols].round(3)
        if len(score_cols) > 1:
            score_cols = ['count'] + score_cols
        sorted_cols = ['target_id', 'target_chromo', 'target_window_start', 'target_window_end',
                       'target_strand', 'sno_id', 'sno_window_start',
                       'sno_window_end'] + score_cols + ['target_seq', 'sno_seq']
        df = df[sorted_cols]
        df.to_csv(tempfile, index=False, mode=mode, header=header)
    os.rename(tempfile, outfile)


def get_target_seq(df_gtf, chromo_dir, target_ids, verbose):
    df_target = df_gtf[df_gtf.feature_id.isin(target_ids)].copy(deep=True)
    df_target = df_target[['seqname', 'start', 'end', 'strand', 'feature_id']]
    df_target = get_seq(df_target, chromo_dir, verbose)
    df_target['seq'] = df_target['seq'].str.upper().str.replace('T', 'U')
    df_target['feature_id'] =  df_target['feature_id'] + '_' + df_target.seqname.astype(str) + '_' + df_target.strand
    df_target = df_target.set_index('feature_id')
    target_dict = df_target.to_dict(orient='index')
    return target_dict


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("sno_fasta", help="fasta file containing snoRNA sequences")
    parser.add_argument("target_ids",
                        help="txt file containing target identifiers (gene_id, transcript_id or exon_id), "
                             "ids should match the ones in the gtf file")
    parser.add_argument("gtf", help="Annotation file in .gtf format. Preferably an annotation of the whole "
                                    "genome or whole chromosomes of specified targets")
    parser.add_argument("chromo_fasta_dir", help="Directory containing fasta files of individual chromosome")
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
    parser.add_argument("--verbose", help='Print the steps', action='store_true')

    if sys.argv[1] in ['-v', '--version']:
        print('snoGloBe version :', __version__)
        exit(0)

    args = parser.parse_args()

    sno_fasta = args.sno_fasta
    target_ids = args.target_ids
    gtf = args.gtf
    chromo_dir  = args.chromo_fasta_dir
    nb_threads = args.nb_threads
    step = args.stepsize
    chunksize = args.chunksize
    threshold = args.threshold
    outfile = os.path.abspath(args.output)
    outpath = os.path.dirname(outfile)
    merge_conswindows = args.merge
    nb_windows = args.nb_windows
    add_seq = args.seq
    verbose = args.verbose

    tests.check_dependencies()
    tests.check_gtf(gtf)
    tests.check_output(outfile)

    if nb_windows != 1 and merge_conswindows is False:
        print('-w/--nb_windows must be used with -m/--merge option. -w/--windows %d will be ignored' % nb_windows,
              file=sys.stderr)

    if threshold < 0 or threshold >1:
        print("Wrong choice of -t/--threshold. Value must be between 0 and 1, chosen value was: ", threshold,
              file=sys.stderr)
        exit(1)

    # define temp file names
    sno_file = os.path.join(outpath, os.path.basename(sno_fasta).replace('.fa', '.input.csv'))
    bedfile = os.path.join(outpath, os.path.basename(gtf).replace('.gtf', '.features.bed'))
    target_file = os.path.join(outpath, os.path.basename(target_ids).replace('.txt', '.input.csv'))

    df_gtf = read_gtf(gtf, verbose)

    with open(target_ids) as tf:
        target_list =  tf.readlines()
        target_list = [t.strip() for t in target_list]
    tests.check_target_ids(df_gtf, target_list)
    sno_windows(sno_fasta, sno_file, verbose)
    if verbose:
        print('Preparing target windows')
    make_bed(df_gtf, bedfile)
    target_dict = get_target_seq(df_gtf, chromo_dir, target_list, verbose)
    target_windows(target_dict, df_gtf, step, target_file, bedfile)
    del df_gtf
    if verbose:
        print('Starting predictions')
    make_pred(step, sno_file, target_file, chunksize, nb_threads, outfile, threshold, merge_conswindows)

    # remove temp files
    os.remove(sno_file)
    os.remove(target_file)
    os.remove(bedfile)

    # merge consecutive windows
    if merge_conswindows is True:
        if verbose:
            print('Merging consecutive windows')
        cons_windows(outfile, nb_windows, chunksize, step, nb_threads)

    if add_seq is True:
        if verbose:
            print('Adding interaction sequences')
        target_dict = {k.rsplit('_', 3)[0]: v for k,v in target_dict.items()}
        interaction_sequence(outfile, sno_fasta, target_dict, chunksize)


if __name__ == '__main__':
    main()
