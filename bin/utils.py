import requests
import pandas as pd
import sys


def fetch_from_ensembl(target_ids):
    server = "https://rest.ensembl.org"
    ext_gtf = "/overlap/id/%s?feature=gene&feature=exon&feature=CDS&feature=transcript"
    ext_seq = "/sequence/id/%s?"

    gtf = []
    d_target = {}
    for t in target_ids:
        r_gtf = requests.get(server + ext_gtf % t, headers={ "Content-Type" : "application/json"})
        r_seq = requests.get(server + ext_seq % t, headers={"Content-Type": "text/plain"})
        if not r_gtf.ok:
            r_gtf.raise_for_status()
            sys.exit(1)
        if not r_seq.ok:
            r_seq.raise_for_status()
            sys.exit(1)
        gtf.extend(r_gtf.json())
        d_target[t] = {'seq': r_seq.text.upper().replace('T', 'U')}

    df_gtf = pd.DataFrame(gtf)

    # build gtf
    df_gtf.loc[df_gtf.strand == -1, 'strand'] = '-'
    df_gtf.loc[df_gtf.strand == 1, 'strand'] = '+'
    df_gtf = df_gtf.rename(columns={'seq_region_name': 'seqname', 'feature_type': 'feature'})
    df_gtf.loc[df_gtf.feature == 'exon', 'transcript_id'] = df_gtf.loc[df_gtf.feature == 'exon', 'Parent']
    df_gtf.loc[df_gtf.feature == 'cds', 'transcript_id'] = df_gtf.loc[df_gtf.feature == 'cds', 'Parent']
    df_gtf.loc[df_gtf.feature == 'transcript', 'gene_id'] = df_gtf.loc[df_gtf.feature == 'transcript', 'Parent']

    df_utr = df_gtf[df_gtf.feature == 'exon'].copy(deep=True)

    if not df_utr.empty:
        d_cds = df_gtf[df_gtf.feature == 'cds'].groupby(['Parent'])[['start', 'end']].agg([min, max]).to_dict()

        df_utr['cds_start'] = df_utr['Parent'].map(d_cds[('start', 'min')])
        df_utr['cds_end'] = df_utr['Parent'].map(d_cds[('end', 'max')])
        df_utr = df_utr.dropna(subset=['cds_start', 'cds_end'])

        df_utr.loc[df_utr.cds_start >= df_utr.start, 'end'] = \
            df_utr.loc[df_utr.cds_start >= df_utr.start][['end', 'cds_start']].min(axis=1)
        df_utr.loc[(df_utr.cds_start >= df_utr.start) & (df_utr.strand == '+'), 'feature'] = 'five_prime_utr'
        df_utr.loc[(df_utr.cds_start >= df_utr.start) & (df_utr.strand == '-'), 'feature'] = 'three_prime_utr'

        df_utr.loc[df_utr.cds_end <= df_utr.end, 'start'] = \
            df_utr.loc[df_utr.cds_end <= df_utr.end][['start', 'cds_end']].max(axis=1)
        df_utr.loc[(df_utr.cds_end <= df_utr.end) & (df_utr.strand == '-'), 'feature'] = 'five_prime_utr'
        df_utr.loc[(df_utr.cds_end <= df_utr.end) & (df_utr.strand == '+'), 'feature'] = 'three_prime_utr'

        df_utr = df_utr[df_utr.feature.str.contains('_prime_utr')]
        df_utr['id'] = df_utr['id'] + '.' + df_utr['feature']

        df_gtf = df_gtf.append([df_utr])
        df_gtf = df_gtf.rename(columns={'id': 'feature_id'})
        df_gtf['score'] = '.'

    # exon and transcripts should inherit gene_biotype (and gene_id?)
    df_gtf.loc[df_gtf.feature == 'gene', 'gene_biotype'] = df_gtf.loc[df_gtf.feature == 'gene', 'biotype']
    df_gtf.loc[df_gtf.feature == 'gene', 'gene_id'] = df_gtf.loc[df_gtf.feature == 'gene', 'feature_id']
    df_gtf.loc[df_gtf.feature == 'transcript', 'transcript_id'] = df_gtf.loc[df_gtf.feature == 'transcript', 'feature_id']
    df_gtf.loc[df_gtf.feature == 'exon', 'exon_id'] = df_gtf.loc[df_gtf.feature == 'exon', 'feature_id']
    d_gene = df_gtf[df_gtf.feature == 'gene'][['gene_id', 'gene_biotype']].set_index('gene_id').to_dict()
    df_gtf.loc[df_gtf.feature == 'transcript', 'gene_biotype'] = df_gtf.loc[df_gtf.feature == 'transcript', 'Parent'].map(d_gene['gene_biotype'])
    d_trx = df_gtf[df_gtf.feature == 'transcript'][['transcript_id', 'gene_id', 'gene_biotype']].set_index('transcript_id').to_dict()
    df_gtf.loc[df_gtf.feature == 'exon', 'gene_biotype'] = df_gtf.loc[df_gtf.feature == 'exon', 'Parent'].map(d_trx['gene_biotype'])
    df_gtf.loc[df_gtf.feature == 'exon', 'gene_id'] = df_gtf.loc[df_gtf.feature == 'exon', 'Parent'].map(d_trx['gene_id'])

    df_gtf[['start', 'end']] = df_gtf[['start', 'end']].astype(int)

    # complete d_target
    for k in target_ids:
        k_start = df_gtf[df_gtf['feature_id'] == k].start.values[0]
        k_strand = df_gtf[df_gtf['feature_id'] == k].strand.values[0]
        k_seqname = df_gtf[df_gtf['feature_id'] == k].seqname.values[0]

        d_target[k]['strand'] = k_strand
        d_target[k]['start'] = k_start
        new_key = str(k_seqname) + '_' + str(k_start) + '_' + k_strand
        d_target[new_key] = d_target.pop(k)

    return df_gtf, d_target


def fetch_from_refseq(target_ids):
    from Bio import Entrez
    # test ids 'NR_046018' 'NR_003693.1' 'NM_001354809'
    for t in target_ids:
        handle = Entrez.esearch(db="nucleotide", term=t)
        rec1 = Entrez.read(handle)
        id = rec1['IdList'][0]
        f = Entrez.efetch(db="nucleotide", id=id, retmode='xml')
        record = Entrez.read(f)


def read_fasta_from_str(fasta_string, seq_type):
    d_fasta = {}
    if type(fasta_string) == str:
        fasta_string = fasta_string.split()

    id = ''
    for line in fasta_string:
        try:
            line = line.decode('utf-8')
        except AttributeError:
            pass
        line = line.strip()
        if line[0] == '>':
            if id != '':
                d_fasta[id]['seq'] = seq
            id = line[1:]
            if seq_type == 'target':
                id += '_0_+'
            d_fasta[id] = {'seq': '', 'start': 0, 'strand': '+'}
            seq = ''
        else:
            seq += line.upper().replace('T', 'U')
    d_fasta[id]['seq'] = seq
    return d_fasta