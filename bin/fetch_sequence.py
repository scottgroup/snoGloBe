import os
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import gzip
import pandas as pd
import gc


def find_prefix(filelist):
    pref = ''
    for i in range(0, len(max(filelist, key=len))):
        pref += filelist[0][i]
        for f in filelist:
            if f[:i+1] != pref:
                return pref[:-1]
    return ''


def find_suffix(filelist):
    suf = ''
    for i in range(1, len(max(filelist, key=len)) + 1):
        suf = filelist[0][-i] + suf
        for f in filelist:
            if f[-i:] != suf:
                return suf[1:]
    return ''


def get_chr_dict(chr_file, chromo, max_pos):
    chr_seq_dict = {}
    # Open an optionally gzipped file
    fn_open = gzip.open if chr_file.endswith('.gz') else open

    file_string = str()
    with fn_open(chr_file) as f:
        n = 0
        for line in f:
            try:
                line = line.decode('utf-8')
            except AttributeError:
                pass
            if line.startswith('>'):
                continue
            else:
                file_string += line.strip()
                n += len(line.strip())
                if n > max_pos:
                    break
    chr_seq_dict[chromo] = file_string
    return chr_seq_dict


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
            sequence = str(Seq(sequence, generic_dna).reverse_complement())
        row['seq'] = sequence
        return row

    dataf['seqname'] = dataf['seqname'].map(str)
    chr_list = list(dataf['seqname'].unique())

    df_out = pd.DataFrame()

    listdir = os.listdir(chr_directory)
    pref = find_prefix(listdir)
    suf = find_suffix(listdir)

    for chr_file in listdir:
        chromo = chr_file.replace(pref, '').replace(suf, '')
        if len(listdir) == 1:
            chromo = chr_file.split('.fa')[0]
        if chromo not in chr_list:
            if 'chr' + chromo in chr_list:
                chromo = 'chr' + chromo
            else:
                continue
        chr_file = os.path.join(chr_directory, chr_file)
        dataf_temp = dataf[dataf.seqname == chromo].copy(deep=True)
        max_pos = dataf_temp.end.max()
        chr_seq_dict = get_chr_dict(chr_file, chromo, max_pos)
        dataf_temp = dataf_temp.apply(fetch_seqs, axis=1)
        df_out = df_out.append(dataf_temp)
        del chr_seq_dict

    gc.collect()
    return df_out
