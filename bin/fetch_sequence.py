import os
from Bio.Seq import Seq
import gzip
import pandas as pd
import gc
import sys


def get_sequence_from_coordinates(dataf, chr_directory, verbose):
    """
    :param dataf: dataframe with columns "seqname", "start", "end" and "strand"
    :param chr_directory : directory holding individual fasta files for each chromosomes.
    :return: same dataframe, but now with a column named "sequence"
    """
    if verbose:
        print('Reading chromosome fasta')

    def fetch_seqs(row):
        sequence = chr_seq_dict[row.seqname][int(row.start - 1):int(row.end)]
        if row.strand == '-':
            sequence = str(Seq(sequence).reverse_complement())
        row['seq'] = sequence
        return row

    dataf['seqname'] = dataf['seqname'].map(str)
    chr_list = list(dataf['seqname'].unique())

    df_out = pd.DataFrame(columns=['seqname', 'start', 'end', 'strand', 'seq'])

    listdir = os.listdir(chr_directory)
    listdir = [f for f in listdir if '.fa' in f]

    for chr_file in listdir:
        chr_file = os.path.join(chr_directory, chr_file)
        chr_seq_dict = {}
        # Open an optionally gzipped file
        fn_open = gzip.open if chr_file.endswith('.gz') else open

        chromo = ''
        with fn_open(chr_file) as f:
            for line in f:
                try:
                    line = line.decode('utf-8')
                except AttributeError:
                    pass
                if line.startswith('>'):
                    if chromo in chr_seq_dict.keys():
                        chr_seq_dict[chromo] = filestring
                    n = 0
                    chromo = line[1:].split()[0]
                    if chromo in chr_list:
                        filestring = str()
                        chr_seq_dict[chromo] = filestring
                        max_pos = dataf[dataf.seqname == chromo].end.max()
                else:
                    if chromo in chr_seq_dict.keys() and n <= max_pos:
                        line = line.strip()
                        filestring += line
                        n += len(line)
        if chromo in chr_seq_dict.keys():
            chr_seq_dict[chromo] = filestring

        dataf_temp = dataf[dataf.seqname.isin(chr_seq_dict.keys())].copy(deep=True)
        dataf_temp = dataf_temp.apply(fetch_seqs, axis=1)
        if not dataf_temp.empty:
            df_out = df_out.append(dataf_temp)
        del chr_seq_dict

        chromo_out = df_out.seqname.unique().tolist()
        missing_chromo = set(chr_list) - set(chromo_out)
        if len(missing_chromo) == 0:
            break

    if len(missing_chromo) > 0:
        print('Error! Missing fasta entry for chromosome(s): %s. Exiting.'% ', '.join(missing_chromo), file=sys.stderr)
        sys.exit(1)

    gc.collect()
    return df_out
