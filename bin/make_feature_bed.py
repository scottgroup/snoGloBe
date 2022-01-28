#!/usr/bin/env python3

import pandas as pd
import subprocess
import os
import sys


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
    df = df.replace({'gene_biotype': biotype_dict})
    return df


def intron_bed(df_gtf):
    df_gtf = df_gtf[df_gtf.feature == 'exon']
    df_gtf = df_gtf.sort_values(['seqname', 'gene_id', 'transcript_id', 'start', 'end'])
    df_gtf = df_gtf[['seqname', 'start', 'end', 'strand', 'transcript_id', 'score']]
    df_gtf['next_start'] = df_gtf.start.shift(-1)
    df_gtf['next_trx'] = df_gtf.transcript_id.shift(-1)
    df_introns = df_gtf[df_gtf.transcript_id == df_gtf.next_trx].copy(deep=True)
    df_introns['intron_start'] = df_introns.end.astype(int) + 1
    df_introns['intron_end'] = df_introns.next_start.astype(int) - 1
    df_introns['feature'] = 'intron'
    return df_introns[['seqname', 'intron_start', 'intron_end', 'feature', 'score', 'strand']]


def make_bed(df_gtf, bedfile, threads, calc_intron=True):
    feature_cols = ['exon', 'five_prime_utr', 'three_prime_utr', 'intron']
    biotype_cols = ['Mt_rRNA', 'Mt_tRNA', 'TEC', 'TR_V_gene',  'lncRNA', 'miRNA', 'misc_RNA',
                    'protein_coding', 'pseudogene', 'rRNA', 'ribozyme', 'sRNA', 'scRNA',
                    'scaRNA', 'snRNA', 'snoRNA', 'tRNA', 'vaultRNA']
    df_biotype = df_gtf[df_gtf.feature == 'gene'][['seqname', 'start', 'end', 'gene_biotype', 'score', 'strand']]
    df_biotype = group_bitoypes(df_biotype)
    df_biotype = df_biotype.drop_duplicates()
    df_biotype = df_biotype.rename(columns={'gene_biotype': 'feature'})
    warn_biotype =  df_biotype[~df_biotype.feature.isin(biotype_cols)].feature.unique().tolist()
    if len(warn_biotype) > 0:
        print('Warning: The following features won\'t be taking into account for the prediction. '
              'Please refer to the manual for a list of accepted features.', file=sys.stderr)
        print(warn_biotype, file=sys.stderr)
    df_biotype = df_biotype[df_biotype.feature.isin(biotype_cols)]

    df_gtf = df_gtf[df_gtf.feature.isin(feature_cols)]
    if calc_intron:
        df_intron = intron_bed(df_gtf)
        df_intron = df_intron.rename(columns={'intron_start': 'start', 'intron_end': 'end'})
    else:
        df_intron = pd.DataFrame()

    df_bed = df_gtf[df_gtf.feature.isin(feature_cols)][['seqname', 'start', 'end', 'feature', 'score', 'strand']]
    df_bed = pd.concat([df_bed, df_intron, df_biotype], sort=False)
    df_bed['end'] += 1  # to bypass the fact that if only the last nucleotide overlaps, the overlap is not reported
    df_bed = df_bed.sort_values(['seqname', 'start', 'end'])
    df_bed[['start', 'end']] = df_bed[['start', 'end']].astype(int)
    df_bed = df_bed.drop_duplicates()

    cols = feature_cols + biotype_cols
    cols.sort()
    for strand in ['+', '-']:
        filelist = []
        for ft in cols:
            tmp_file = bedfile + '.' + ft + strand + '.temp.bed'
            outfile = bedfile + '.' + ft + strand + '.bed'
            filelist.append(outfile)
            df_temp = df_bed[(df_bed.feature == ft) & (df_bed['strand'] == strand)]
            df_temp.to_csv(tmp_file , sep='\t', header=False, index=False)
            if not df_temp.empty:
                bedtools_cmd = ('bedtools', 'merge', '-s', '-c', '4,5,6', '-o', 'distinct', '-i', tmp_file)
                res = subprocess.Popen(bedtools_cmd, stdout=open(outfile, mode='a'))
                res.wait()
                os.remove(tmp_file)
            else:
                os.rename(tmp_file, outfile)

        # split by strand
        multiint_cmd = ('bedtools', 'multiinter', '-i') + tuple(filelist) #+ ('-header', '-names') + tuple(os.path.basename(i).rsplit('.',2)[1] for i in filelist)
        awk_cmd = ("awk", "{print $" + '"\t"$'.join([str(i) for i in range(1, 4)]) + '"\t.\t.\t' + strand
                   + '\t"$' + '"\t"$'.join([str(i) for i in range(6, len(cols) + 6 + 1)]) + "}")

        res_multi = subprocess.Popen(multiint_cmd, stdout=subprocess.PIPE)
        res_awk = subprocess.Popen(awk_cmd, stdin=res_multi.stdout, stdout=open(bedfile + '.presort', mode='a'))
        res_awk.wait()
        for f in filelist:
            os.remove(f)
    sort_cmd = ('sort', '-k1,1', '-k2,3n', '-k6,6', '--parallel=' + str(threads),  bedfile + '.presort')
    res_sort = subprocess.Popen(sort_cmd, stdout=open(bedfile + '.sorted', mode='a'))
    res_sort.wait()
    os.remove(bedfile + '.presort')
    groupby_cmd = ('bedtools', 'groupby', '-i', bedfile + '.sorted', '-g', '1-6',
                   '-c', ','.join(str(i) for i in range(7, len(cols) + 6 + 1)), '-o', 'max')
    res_grp = subprocess.Popen(groupby_cmd, stdout=open(bedfile, mode='a'))
    res_grp.wait()
    if res_grp.returncode == 0:
        pass
        os.remove(bedfile + '.sorted')
    else:
        exit(1)

