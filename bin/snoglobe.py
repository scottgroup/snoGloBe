#!/usr/bin/env python3

__author__ = "Gabrielle Deschamps-Francoeur"
__email__ = "gabrielle.deschamps-francoeur@usherbrooke.ca"
__version__ = '0.1.1'

import os, sys
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
import argparse
from sliding_windows import sliding_window, read_trx_fasta
from make_feature_bed import make_bed
import pandas as pd
import extract_features as ef
from predict_mp import make_pred
from consecutive_windows import cons_windows
import tests
from fetch_sequence import get_seq_bedtools, get_window_seq
import re
import io
import subprocess
import gc
import random
from utils import fetch_from_ensembl, read_fasta_from_str

pd.set_option('expand_frame_repr', False)


def parse_args():
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
    return args


def read_gtf(gtf_file):
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


def interaction_sequence(outfile, sno_dict, chromo_dict, chunksize, cols, temp_dir):
    tempfile1 = temp_dir + '.tseq'
    get_window_seq(chromo_dict, outfile, tempfile1)
    for i, df in enumerate(pd.read_csv(tempfile1, chunksize=chunksize, names=cols, sep='\t')):
        if 'add_seq' in df.columns:
            df = df.drop(['add_seq'], axis=1)
        if i == 0:
            mode = 'w'
        else:
            mode = 'a'
        if df.empty is False:
            df['sno_seq'] = df.apply(lambda row: extract_seq(row, sno_dict, 'sno'), axis=1)
            df['target_seq'] = df['target_seq'].str.upper().str.replace('T', 'U')
        else:
            df['sno_seq'] = None
            df['target_seq'] = None
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
        df.to_csv(outfile, index=False, mode=mode, header=False, sep='\t')
    os.remove(tempfile1)


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


def get_target_seq(df_gtf, chromo_dir, target_ids, verbose, outpath, step):
    df_target = df_gtf[df_gtf.feature_id.isin(target_ids)].copy(deep=True)
    df_target = df_target[['seqname', 'start', 'end', 'strand', 'feature_id']]
    chromo_dict = get_seq_bedtools(df_target, chromo_dir, verbose, outpath, step)
    return chromo_dict


def run_snoglobe(sno_dict, target_list, chromo_dir, df_gtf, nb_threads, step, chunksize, threshold,
                 outfile, temp_dir, nb_windows, merge_conswindows=True, add_seq=True, verbose=False, calc_intron=True):

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
    chromo_dict = get_target_seq(df_gtf, chromo_dir, target_list, verbose, temp_dir, step)
    make_bed(df_gtf, bedfile, nb_threads, calc_intron)
    ef.extract_features(temp_dir, bedfile, target_file)

    del df_gtf
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
        cons_windows(outfile, nb_windows, chunksize, step, nb_threads, add_seq)
        cols = ['target_chromo', 'target_window_start', 'target_window_end', 'target_id', 'count', 'target_strand',
                'sno_id', 'sno_window_start', 'sno_window_end', 'mean_score', 'min_score', 'max_score']

    if add_seq:
        if verbose:
            print('Adding interaction sequences')
        if merge_conswindows:
            cols.append('add_seq')
        cols.extend(['target_seq'])
        interaction_sequence(outfile, sno_dict, chromo_dict, chunksize, cols, temp_dir)
        if merge_conswindows:
            cols.remove('add_seq')
        cols.extend(['sno_seq'])

    temp_chromo_files = [temp_dir + '.chromo.fa', temp_dir + '.chromo.fa.fai']
    for temp_f in temp_chromo_files:
        if os.path.exists(temp_f):
            os.remove(temp_f)

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
        res2.wait()
        res1.wait()
        res0.wait()

        if res0.returncode == 0 and res1.returncode == 0 and res2.returncode == 0 and res3.returncode == 0:
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
    args = parse_args()

    sno_fasta = args.sno_fasta
    target_ids = args.target_ids
    full_gtf = args.gtf
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
    tests.check_gtf(full_gtf)
    tests.check_output(outfile)

    if not os.path.isfile(sno_fasta):
        print("snoRNA fasta file could not be opened :", sno_fasta, '. Exiting.', file=sys.stderr)
        exit(1)
    if not os.path.isfile(target_ids):
        print("Target id file could not be opened :", target_ids, '. Exiting.', file=sys.stderr)
        exit(1)

    if nb_windows != 1 and merge_conswindows is False:
        print('-w/--nb_windows must be used with -m/--merge option. -w/--windows %d will be ignored' % nb_windows,
              file=sys.stderr)

    if threshold < 0 or threshold >1:
        print("Wrong choice of -t/--threshold. Value must be between 0 and 1, chosen value was: ", threshold,
              file=sys.stderr)
        exit(1)

    comb_run = 'temp_' + str(random.randint(0, 100000))
    temp_dir = os.path.join(os.path.dirname(outfile), comb_run, os.path.basename(outfile))
    os.makedirs(os.path.dirname(temp_dir), exist_ok=True)

    if verbose:
        print('Reading gtf')
    # use smallest portion of gtf possible
    grep_cmd = ('grep', '-f', target_ids, full_gtf)
    res0 = subprocess.Popen(grep_cmd, stdout=open(temp_dir + '.target.gtf', 'a'))
    res0.wait()
    if res0.returncode != 0:
        sys.exit(1)

    # bedtools intersect target_bed gtf > small.gtf
    gtf = temp_dir + '.small.gtf'
    bedtools_cmd = ('bedtools', 'intersect', '-s', '-u', '-a', full_gtf, '-b', temp_dir + '.target.gtf')
    res1 = subprocess.Popen(bedtools_cmd, stdout=open(gtf, 'a'))
    res1.wait()
    if res1.returncode == 0:
        os.remove(temp_dir + '.target.gtf')
    else:
        sys.exit(1)

    # read small gtf
    df_gtf = read_gtf(gtf)
    os.remove(gtf)

    with open(target_ids) as tf:
        target_list =  tf.readlines()
        target_list = [t.strip() for t in target_list]
    tests.check_target_ids(df_gtf, target_list)


    sno_dict = read_trx_fasta(sno_fasta)

    run_snoglobe(sno_dict, target_list, chromo_dir, df_gtf, nb_threads, step, chunksize, threshold,
                 outfile, temp_dir, nb_windows, merge_conswindows=merge_conswindows, add_seq=add_seq,
                 verbose=verbose, calc_intron=True)


def prep_snoglobe_web(step, threshold, nb_windows, nb_threads, chunksize, outfile, sno_fasta='', target_ids='',
                      target_fasta='', target_info='', target_biotype=''):

    temp_dir = os.path.join(outfile.rsplit('.',1)[0], os.path.basename(outfile))
    os.makedirs(os.path.dirname(temp_dir))

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

    with open(temp_dir + '.chromo.fa', 'w') as w:
        for k in target_dict.keys():
            w.write('>' + k.rsplit('_',2)[0] + '\n' + target_dict[k]['seq'] + '\n')
    del target_dict

    sno_dict = read_fasta_from_str(sno_fasta, 'sno')

    run_snoglobe(sno_dict, target_list, os.path.dirname(temp_dir), df_gtf, nb_threads, step,
                 chunksize, threshold, outfile, temp_dir, nb_windows, calc_intron=calc_intron, verbose=True)


if __name__ == '__main__':
    prep_snoglobe()

