#!/usr/bin/env

import pickle
import pandas as pd
import multiprocessing as mp
import numpy as np
import gc
import os
from extract_features import set_columns
from sliding_windows import seq_onehotmatrix


def predict(args):
    X_sno, X_t, mdl, step, threshold = args
    X_sno['decoy'] = 0
    X_t['decoy'] = 0
    X = pd.merge(X_sno, X_t, on='decoy')
    X = X.drop(['decoy'], axis=1)
    X = X.set_index(X.snoid + ',' + X.wid)
    X = X.drop(['snoid', 'wid'], axis=1)
    X = X[X.columns.sort_values()]
    X['rel_pos'] = X['rel_pos'].round(3)
    res = mdl.predict_proba(X)
    df_proba = pd.DataFrame(res, index=X.index, columns=['prob_neg', 'prob_pos'])
    df_proba['snoid'] = df_proba.index.str.split(',', expand=True).get_level_values(0)
    df_proba['target_id'] = df_proba.index.str.split(',', expand=True).get_level_values(1)
    df_proba[['snoid', 'sno_window_start']] = df_proba['snoid'].str.rsplit('_', 1, expand=True)
    df_proba['sno_window_end'] = df_proba['sno_window_start'].astype(int) + 13
    df_proba[['target_chromo', 'target_strand', 'wstart']] = df_proba['target_id'].str.rsplit('_', 2, expand=True)
    df_proba['target_id'] = df_proba['target_id'].str.rsplit('_', 1).str[0]
    df_proba['wend'] = df_proba.wstart.astype(int) + 13
    df_proba = df_proba[['target_chromo', 'wstart', 'wend', 'target_id', 'prob_pos', 'target_strand',
                         'snoid', 'sno_window_start', 'sno_window_end']]
    df_proba = df_proba[df_proba.prob_pos > threshold]
    df_proba.prob_pos = df_proba.prob_pos.round(3)
    return df_proba


def make_pred(step, snofile, targetfile, chunksize, nb_threads, outfile, threshold):
    snocols = [
        'snoid', '0_A', '0_C', '0_G', '0_U', '10_A', '10_C', '10_G', '10_U', '11_A', '11_C', '11_G', '11_U',
        '12_A', '12_C', '12_G', '12_U', '1_A', '1_C', '1_G', '1_U', '2_A', '2_C', '2_G', '2_U', '3_A', '3_C',
        '3_G', '3_U', '4_A', '4_C', '4_G', '4_U', '5_A', '5_C', '5_G', '5_U', '6_A', '6_C', '6_G', '6_U',
        '7_A', '7_C', '7_G', '7_U', '8_A', '8_C', '8_G', '8_U', '9_A', '9_C', '9_G', '9_U', 'rel_pos'
    ]
    targetcols = [
        'wid', 'seq', 'Mt_rRNA', 'Mt_tRNA', 'TEC', 'TR_V_gene', 'exon', 'five_prime_utr',
        'intron', 'lncRNA', 'miRNA', 'misc_RNA', 'protein_coding', 'pseudogene', 'rRNA', 'ribozyme', 'sRNA',
        'scRNA', 'scaRNA', 'snRNA', 'snoRNA', 'tRNA', 'three_prime_utr', 'vaultRNA'
    ]
    cat_tcol = targetcols.copy()
    cat_snocol = snocols.copy()
    cat_snocol.remove('rel_pos')
    cat_snocol.remove('snoid')
    cat_tcol.remove('wid')
    cat_tcol.remove('seq')

    snodtypes = {i: "bool" for i in cat_snocol}
    tdtypes = {i: "bool" for i in cat_tcol}

    nb_snowindows = 0
    with open(snofile, 'r') as f:
        for line in f:
            nb_snowindows += 1

    if nb_snowindows < 500:
        if chunksize < nb_snowindows:
            chunksize_sno = chunksize
        else:
            chunksize_sno = nb_snowindows
        chunksize_t = int(chunksize / chunksize_sno)
    else:
        chunksize_sno = 500
        chunksize_t = int(chunksize / 500)

    mdl_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'snoGloBe.sav')
    mdl = pickle.load(open(mdl_path, 'rb'))
    mode = 'w'
    for X_t in pd.read_csv(targetfile, names=targetcols, dtype=tdtypes, chunksize=chunksize_t, sep='\t'):
        X_t['seq'] = X_t['seq'].str.upper().str.replace('T', 'U')
        X_t = seq_onehotmatrix(X_t)
        X_t = set_columns(X_t)
        tseq_cols = [c for c in  X_t.columns if c[0].isdigit()]
        X_t[tseq_cols] = X_t[tseq_cols].astype(bool)
        for X_sno in pd.read_csv(snofile, dtype=snodtypes, chunksize=chunksize_sno):
            for i in set(snocols) - set(X_sno.columns):
                X_sno[i] = False
            X_sno = X_sno[snocols]
            pool = mp.Pool(processes=nb_threads)
            res = pool.map(predict,
                           [[X_sno, X_part, mdl, step, threshold] for X_part in np.array_split(X_t, nb_threads)])
            pool.close()
            df_proba = pd.concat(list(res))
            df_proba.to_csv(outfile, header=False, index=False, mode=mode, sep='\t')
            mode = 'a'
            del df_proba
            gc.collect()
