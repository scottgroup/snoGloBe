#!/usr/bin/env python3

import pandas as pd


def intron_bed(df_gtf):
    df_gtf = df_gtf[df_gtf.feature == 'exon']
    df_gtf = df_gtf[['seqname', 'start', 'end', 'strand', 'transcript_id', 'score']]
    df_gtf['next_start'] = df_gtf.start.shift(-1)
    df_gtf['next_trx'] = df_gtf.transcript_id.shift(-1)
    df_introns = df_gtf[df_gtf.transcript_id == df_gtf.next_trx].copy(deep=True)
    df_introns.loc[df_introns.strand == '+', 'intron_start'] = df_introns.end
    df_introns.loc[df_introns.strand == '+', 'intron_end'] = df_introns.next_start
    df_introns.loc[df_introns.strand == '-', 'intron_end'] = df_introns.end
    df_introns.loc[df_introns.strand == '-', 'intron_start'] = df_introns.next_start
    df_introns[['intron_end']] = df_introns[['intron_end']].astype(int)
    df_introns['intron_end'] -= 1
    df_introns['intron_start'] += 1
    df_introns['feature'] = 'intron'
    return df_introns[['seqname', 'intron_start', 'intron_end', 'feature', 'score', 'strand']]


def make_bed(df_gtf, bedfile):
    df_bed = df_gtf[df_gtf.feature.isin(
        ['five_prime_utr', 'three_prime_utr', 'exon'])][['seqname', 'start', 'end', 'feature', 'score', 'strand']]
    df_intron = intron_bed(df_gtf)
    df_intron = df_intron.rename(columns={'intron_start': 'start', 'intron_end': 'end'})
    df_biotype = df_gtf[df_gtf.feature == 'gene'][['seqname', 'start', 'end', 'gene_biotype', 'score', 'strand']]
    df_biotype = df_biotype.rename(columns={'gene_biotype': 'feature'})
    df_bed = pd.concat([df_bed, df_intron, df_biotype], sort=False)
    df_bed['end'] += 1  # to bypass the fact that if only the last nucleotide overlaps, the overlap is not reported
    df_bed = df_bed.sort_values(['seqname', 'start', 'end'])
    df_bed[['start', 'end']] = df_bed[['start', 'end']].astype(int)
    df_bed.to_csv(bedfile, sep='\t', header=False, index=False)
