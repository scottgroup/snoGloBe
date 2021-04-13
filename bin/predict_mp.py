#!/usr/bin/env

import pickle
import pandas as pd
import multiprocessing as mp
import numpy as np
import gc
import os


def predict(args):
    X_sno, X_t, mdl, step, threshold = args
    X_sno['decoy'] = 0
    X_t['decoy'] = 0
    X = pd.merge(X_sno, X_t, on='decoy')
    del X_t
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
    df_proba[['target_chromo', 'start', 'target_strand', 'wstart']] = df_proba['target_id'].str.rsplit('_', 3,
                                                                                                       expand=True)
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
        '13_A', '13_C', '13_G', '13_U', '23_A', '23_C', '23_G', '23_U', '24_A', '24_C', '24_G', '24_U',
        '25_A', '25_C', '25_G', '25_U', '14_A', '14_C', '14_G', '14_U', '15_A', '15_C', '15_G', '15_U',
        '16_A', '16_C', '16_G', '16_U', '17_A', '17_C', '17_G', '17_U', '18_A', '18_C', '18_G', '18_U',
        '19_A', '19_C', '19_G', '19_U', '20_A', '20_C', '20_G', '20_U', '21_A', '21_C', '21_G', '21_U',
        '22_A', '22_C', '22_G', '22_U', 'Mt_rRNA', 'Mt_tRNA', 'TEC', 'TR_V_gene', 'exon', 'five_prime_utr',
        'intron', 'lncRNA', 'miRNA', 'misc_RNA', 'protein_coding', 'pseudogene', 'rRNA', 'ribozyme', 'sRNA',
        'scRNA', 'scaRNA', 'snRNA', 'snoRNA', 'tRNA', 'three_prime_utr', 'vaultRNA', 'wid'
    ]
    cat_tcol = targetcols.copy()
    cat_snocol = snocols.copy()
    cat_snocol.remove('rel_pos')
    cat_snocol.remove('snoid')
    cat_tcol.remove('wid')

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
    for X_sno in pd.read_csv(snofile, dtype=snodtypes, chunksize=chunksize_sno):
        for i in set(snocols) - set(X_sno.columns):
            X_sno[i] = False
        X_sno = X_sno[snocols]
        for X_t in pd.read_csv(targetfile, dtype=tdtypes, chunksize=chunksize_t):
            for i in set(targetcols) - set(X_t.columns):
                X_t[i] = False
            X_t = X_t[targetcols]
            pool = mp.Pool(processes=nb_threads)
            res = pool.map(predict,
                           [[X_sno, X_part, mdl, step, threshold] for X_part in np.array_split(X_t, nb_threads)])
            pool.close()
            df_proba = pd.concat(list(res))
            df_proba.to_csv(outfile, header=False, index=False, mode=mode, sep='\t')
            mode = 'a'
            del df_proba
            gc.collect()
