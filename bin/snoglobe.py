#!/usr/bin/env python3

__author__ = "Gabrielle Deschamps-Francoeur"
__email__ = "gabrielle.deschamps-francoeur@usherbrooke.ca"
__version__ = '0.1.1'

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
import io
import subprocess
import pickle
import gc
import random
from utils import fetch_from_ensembl, read_fasta_from_str


def read_gtf(gtf_file, verbose=False):
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


def sno_windows(sno_dict, sno_file, verbose):
    # read sno fasta file and split in sliding windows
    if verbose:
        print('preparing snoRNA windows')
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


def target_windows(target_dict, df_gtf, step, target_file, bedfile, temp_output):
    df_target = sliding_window(target_dict, step)
    df_target['feature_id'] = df_target.index
    df_target[['seqname', 'start', 'strand', 'window_start']] = df_target.feature_id.str.rsplit('_', 3, expand=True)
    df_target['window_end'] = df_target.window_start.astype(int) + 13
    df_target = df_target.drop(['idx', 'length', 'feature_id', 'start'], axis=1)
    df_bed = df_target[['seqname', 'window_start', 'window_end', 'strand']].copy(deep=True)
    df_target = df_target.drop(['seqname', 'strand', 'window_start', 'window_end'], axis=1)
    df_bed['score'] = 0
    df_bed['idx'] = df_bed.index
    df_out = ef.extract_features(df_bed, temp_output, bedfile)
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


def interaction_sequence(outfile, sno_dict, target_dict, chunksize, cols, temp_dir):
    tempfile = temp_dir + '.seq'
    for i, df in enumerate(pd.read_csv(outfile, chunksize=chunksize, names=cols, sep='\t')):
        if i == 0:
            mode = 'w'
        else:
            mode = 'a'
        if df.empty is False:
            df['target_seq'] = df.apply(lambda row: extract_seq(row, target_dict, 'target'), axis=1)
            df['sno_seq'] = df.apply(lambda row: extract_seq(row, sno_dict, 'sno'), axis=1)
        else:
            df['target_seq'] = None
            df['sno_seq'] = None
        score_cols = [c for c in df.columns if 'score' in c]
        df[score_cols] = df[score_cols].round(3)
        if len(score_cols) > 1:
            sorted_cols = ['target_chromo', 'target_window_start', 'target_window_end',
                            'target_id', 'count', 'target_strand', 'sno_id', 'sno_window_start',
                           'sno_window_end'] + score_cols + ['target_seq', 'sno_seq']
        else:
            sorted_cols = ['target_chromo', 'target_window_start', 'target_window_end',
                            'target_id'] + score_cols + ['target_strand', 'sno_id', 'sno_window_start',
                           'sno_window_end', 'target_seq', 'sno_seq']
        df = df[sorted_cols]
        df.to_csv(tempfile, index=False, mode=mode, header=False, sep='\t')
    os.rename(tempfile, outfile)


def merge_intervals(df, strand, outpath):
    temp_file = outpath + '.temp%s.bed' % strand
    # split by strand
    temp_df = df[df.strand == strand][['seqname', 'start', 'end']]
    temp_df = temp_df.sort_values(['seqname', 'start', 'end'])
    temp_df.to_csv(temp_file, index=False, header=False, sep='\t')
    # merge intervals with bedtools
    cmd = ('bedtools',  'merge', '-i', temp_file)
    res = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    df_merged = pd.read_csv(io.StringIO(res.stdout.read().decode()),
                            sep='\t',
                            names=['seqname', 'start', 'end'],
                            dtype={'seqname':str, 'start': int, 'end': int})
    df_merged['strand'] = strand
    os.remove(temp_file)
    return df_merged


def get_target_seq(df_gtf, chromo_dir, target_ids, verbose, outpath):
    df_target = df_gtf[df_gtf.feature_id.isin(target_ids)].copy(deep=True)
    df_target = df_target[['seqname', 'start', 'end', 'strand', 'feature_id']]
    df_fwd = merge_intervals(df_target, '+', outpath)
    df_rev = merge_intervals(df_target, '-', outpath)
    df_merged = pd.concat([df_fwd, df_rev])
    df_merged = get_seq(df_merged, chromo_dir, verbose)
    df_merged['seq'] = df_merged['seq'].str.upper().str.replace('T', 'U')
    df_merged['feature_id'] = df_merged.seqname.astype(str) \
                              + '_' + df_merged.start.astype(str) \
                              + '_' + df_merged.strand
    df_merged = df_merged.set_index('feature_id')
    target_dict = df_merged.to_dict(orient='index')
    return target_dict


def run_snoglobe(sno_dict, target_list, target_dict, df_gtf, nb_threads, step, chunksize, threshold,
                 outfile, temp_dir, nb_windows, merge_conswindows=True, add_seq=True, verbose=False, calc_intron=True):

    outpath = os.path.dirname(outfile)
    # define temp file names

    sno_file = temp_dir + '.snoinput.csv'
    bedfile = temp_dir + '.features.bed'
    target_file = temp_dir + '.tinput.csv'
    target_bed = temp_dir + '.target.bed'

    df_gtf[df_gtf.feature_id.isin(target_list)][
        ['seqname', 'start', 'end', 'feature_id', 'score', 'strand']
    ].to_csv(target_bed, index=False, header=False, sep='\t')

    sno_windows(sno_dict, sno_file, verbose)

    if verbose:
        print('Preparing target windows')
    make_bed(df_gtf, bedfile, calc_intron)
    target_windows(target_dict, df_gtf, step, target_file, bedfile, temp_dir)

    if add_seq:
        target_pickle = temp_dir + '.seq.pkl'
        pickle.dump(target_dict, open(target_pickle, 'wb'))

    del target_dict, df_gtf
    gc.collect()

    if verbose:
        print('Starting predictions')
    make_pred(step, sno_file, target_file, chunksize, nb_threads, outfile, threshold)

    # remove temp files
    os.remove(sno_file)
    os.remove(target_file)
    os.remove(bedfile)

    cols = ['target_chromo', 'target_window_start', 'target_window_end', 'target_id',
            'score', 'target_strand', 'sno_id', 'sno_window_start', 'sno_window_end']

    # merge consecutive windows
    if merge_conswindows:
        if verbose:
            print('Merging consecutive windows')
        cons_windows(outfile, nb_windows, chunksize, step, nb_threads)
        cols = ['target_chromo', 'target_window_start', 'target_window_end', 'target_id', 'count', 'target_strand',
                'sno_id', 'sno_window_start', 'sno_window_end', 'mean_score', 'min_score', 'max_score']

    if add_seq:
        target_dict = pickle.load(open(target_pickle, 'rb'))
        if verbose:
            print('Adding interaction sequences')
        interaction_sequence(outfile, sno_dict, target_dict, chunksize, cols, temp_dir)
        cols.extend(['target_seq', 'sno_seq'])
        os.remove(target_pickle)

    # verify if there are any predicted interactions
    nb_line = 0
    with open(outfile) as f:
        for line in f:
            nb_line += 1
            if nb_line > 0:
                break

    if nb_line > 0:
        # sort bed
        sort_cmd = ('sort', '--parallel', str(nb_threads), '-k', '1,1', '-k', '2,3n', outfile)
        # bedtools intersect to get overlapping target_id
        intersect_cmd = ('bedtools', 'intersect', '-s', '-wa', '-wb', '-a', '-', '-b', target_bed)

        groupby_cmd = (
        'bedtools', 'groupby',
        '-g', ','.join([str(i) for i in range(1, len(cols) + 1)]),
        '-c', str(len(cols) + 4),
        '-o', 'distinct')

        awk_cmd = ("awk", "{print $" + '"\t"$'.join([str(i) for i in range(1, 4)]) + '"\t"$' + str(
            len(cols) + 1) + '"\t"$' + '"\t"$'.join([str(i) for i in range(5, len(cols) + 1)]) + "}")


        final_file = outfile + '.final'
        # create file with header
        print('\t'.join(cols), file=open(final_file, 'w'))

        res0 = subprocess.Popen(sort_cmd, stdout=subprocess.PIPE)
        res1 = subprocess.Popen(intersect_cmd, stdin=res0.stdout, stdout=subprocess.PIPE)
        res2 = subprocess.Popen(groupby_cmd, stdin=res1.stdout, stdout=subprocess.PIPE)
        res3 = subprocess.Popen(awk_cmd, stdin=res2.stdout, stdout=open(final_file, 'a'))
        res3.wait()

        if res3.returncode == 0:
            # clean up temp files
            os.remove(target_bed)
            os.rename(final_file, outfile)
        else:
            sys.exit(1)

    else:
        print('\t'.join(cols), file=open(outfile, 'w'))
        os.remove(target_bed)

    os.rmdir(os.path.dirname(temp_dir))


def prep_snoglobe():
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
    parser.add_argument("-c", "--chunksize", help="Number of combinations to run at once. Default: 2500000",
                        type=int, default=2500000)
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

    df_gtf = read_gtf(gtf, verbose)

    with open(target_ids) as tf:
        target_list =  tf.readlines()
        target_list = [t.strip() for t in target_list]
    tests.check_target_ids(df_gtf, target_list)

    comb_run = 'temp_' + str(random.randint(0, 100000))
    temp_dir = os.path.join(os.path.dirname(outfile), comb_run, os.path.basename(outfile))
    os.makedirs(os.path.dirname(temp_dir))
    sno_dict = read_trx_fasta(sno_fasta)
    target_dict = get_target_seq(df_gtf, chromo_dir, target_list, verbose, temp_dir)

    run_snoglobe(sno_dict, target_list, target_dict, df_gtf, nb_threads, step, chunksize, threshold,
                 outfile, temp_dir, nb_windows, merge_conswindows, add_seq, verbose)


def prep_snoglobe_web(step, threshold, nb_windows, nb_threads, chunksize, outfile, sno_fasta='', target_ids='',
                      target_fasta='', target_info='', target_biotype=''):
    if target_ids != '':
        target_list = target_ids.split()
        df_gtf, target_dict = fetch_from_ensembl(target_list)
        calc_intron = True
    elif target_fasta != '':
        target_dict = read_fasta_from_str(target_fasta, 'target')
        gtf = []
        for k in target_dict.keys():
            k_entry = k.rsplit('_', 2)
            k_len = len(target_dict[k]['seq'])
            for t_loc in target_info:
                gtf.append([k.rsplit('_',2)[0]] + k_entry + [k_len, t_loc, target_biotype, '.'])
        df_gtf = pd.DataFrame(gtf,
                              columns=['feature_id', 'seqname', 'start', 'strand', 'end',
                                       'feature', 'gene_biotype', 'score'])
        df_gtf[['start', 'end']] = df_gtf[['start', 'end']].astype(int)
        target_list = df_gtf.feature_id.unique().tolist()
        df_gtf['gene_id'] = df_gtf['seqname']
        df_gtf['transcript_id'] = df_gtf['seqname'] + '_trx'
        calc_intron = False

    else:
        sys.exit()

    sno_dict = read_fasta_from_str(sno_fasta, 'sno')
    comb_run = 'temp_' + str(random.randint(0, 100000))
    temp_dir = os.path.join(os.path.dirname(outfile), comb_run, os.path.basename(outfile))
    os.makedirs(os.path.dirname(temp_dir))

    run_snoglobe(sno_dict, target_list, target_dict, df_gtf, nb_threads, step, chunksize, threshold,
                      outfile, temp_dir, nb_windows, calc_intron=calc_intron)


if __name__ == '__main__':
    prep_snoglobe()

