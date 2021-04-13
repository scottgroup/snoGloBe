#!/usr/bin/env python3

import pandas as pd
import os
import subprocess
from io import StringIO
import sys


def column_list(df_gtf):
    cols = []
    length = 13
    nts = ['A', 'U', 'C', 'G']
    seq_cols = [str(i) + '_' + str(j) for i in range(0, int(length)) for j in nts]
    cols.extend(seq_cols)
    df_gtf_features = df_gtf[['feature', 'gene_biotype']]
    target_features = ['five_prime_utr', 'three_prime_utr', 'exon', 'intron']
    target_features.extend(df_gtf_features[df_gtf_features.feature == 'gene'].gene_biotype.unique().tolist())
    grouped_biotypes = [
        'IG_C_gene', 'IG_D_gene', 'IG_gene', 'IG_J_gene', 'IGLV_gene', 'IGM_gene', 'IG_V_gene',
        'IGZ_gene', 'nonsense_mediated_decay', 'nontranslating_CDS', 'non_stop_decay',
        'polymorphic_pseudogene', 'TR_C_gene', 'TR_D_gene', 'TR_J_gene', 'disrupted_domain',
        'IG_C_pseudogene', 'IG_J_pseudogene', 'IG_pseudogene', 'IG_V_pseudogene',
        'processed_pseudogene', 'transcribed_processed_pseudogene', 'transcribed_unitary_pseudogene',
        'transcribed_unprocessed_pseudogene', 'translated_processed_pseudogene', 'TR_J_pseudogene',
        'TR_V_pseudogene', 'unitary_pseudogene', 'unprocessed_pseudogene', '3prime_overlapping_ncrna',
        'ambiguous_orf', 'antisense', 'antisense_RNA', 'lincRNA',
        'ncrna_host', 'processed_transcript', 'sense_intronic', 'sense_overlapping',
        'bidirectional_promoter_lncRNA', '3prime_overlapping_ncRNA', 'macro_lncRNA', 'non_coding'
    ]
    target_features = [i for i in target_features if i not in grouped_biotypes]
    target_features.append('lncRNA')
    cols.extend(target_features)
    cols.sort()
    return cols


def group_bitoypes(df):
    protein_coding = ['IG_C_gene', 'IG_D_gene', 'IG_gene', 'IG_J_gene', 'IGLV_gene', 'IGM_gene', 'IG_V_gene',
                      'IGZ_gene', 'nonsense_mediated_decay', 'nontranslating_CDS', 'non_stop_decay',
                      'polymorphic_pseudogene', 'TR_C_gene', 'TR_D_gene', 'TR_J_gene']
    pseudogene = ['disrupted_domain', 'IG_C_pseudogene', 'IG_J_pseudogene', 'IG_pseudogene', 'IG_V_pseudogene',
                  'processed_pseudogene', 'transcribed_processed_pseudogene', 'transcribed_unitary_pseudogene',
                  'transcribed_unprocessed_pseudogene', 'translated_processed_pseudogene', 'TR_J_pseudogene',
                  'TR_V_pseudogene', 'unitary_pseudogene', 'unprocessed_pseudogene']
    lncRNA = ['3prime_overlapping_ncrna', 'ambiguous_orf', 'antisense', 'antisense_RNA', 'lincRNA',
              'ncrna_host', 'processed_transcript', 'sense_intronic', 'sense_overlapping',
              'bidirectional_promoter_lncRNA', '3prime_overlapping_ncRNA', 'macro_lncRNA', 'non_coding']

    biotype_dict = {}
    for i in protein_coding:
        biotype_dict[i] = 'protein_coding'
    for i in pseudogene:
        biotype_dict[i] = 'pseudogene'
    for i in lncRNA:
        biotype_dict[i] = 'lncRNA'
    df = df.replace({'target_feature': biotype_dict})
    return df


def set_columns(df, cols):
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


def extract_features(df, fpath, bedfile):
    df = df[['seqname', 'window_start', 'window_end', 'idx', 'score', 'strand']]
    df = df.sort_values(['seqname', 'window_start', 'window_end'])
    df[['window_start', 'window_end']] = df[['window_start', 'window_end']].astype(int)
    df.to_csv(os.path.join(fpath, '_targets.bed'), sep='\t', index=False, header=False)
    bedtools_cmd = ['bedtools', 'intersect', '-a', os.path.join(fpath, '_targets.bed'),
                    '-b', bedfile, '-s', '-wb']
    result = subprocess.Popen(bedtools_cmd, stdout=subprocess.PIPE)
    b = StringIO(result.communicate()[0].decode('utf-8'))
    df_out = pd.read_csv(b, sep="\t",
                         names=['chr1', 'start1', 'end1', 'id', 'score1', 'strand1',
                                'chr2', 'start2', 'end2', 'target_feature', 'score2', 'strand2'],
                         dtype={'chr1': str, 'chr2': str})
    df_out = group_bitoypes(df_out)
    df_dummy = pd.get_dummies(df_out[['target_feature']])
    df_dummy.columns = [i.replace('target_feature_', '') for i in df_dummy.columns.tolist()]
    df_merged = pd.merge(df_out[['id']], df_dummy, left_index=True, right_index=True)
    df_merged = df_merged.groupby('id').max().reset_index()
    df_merged = df_merged.set_index('id')
    os.remove(os.path.join(fpath, '_targets.bed'))
    return df_merged
