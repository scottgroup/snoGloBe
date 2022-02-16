import os
import sys
import subprocess


def get_chr_file_dict(chr_directory, target_chromo_list):
    listdir = os.listdir(chr_directory)
    listdir = [f for f in listdir if '.fa' in f and not f.endswith('.fai')]
    chromo_file_dict = {}
    for chr_file in listdir:
        chr_file = os.path.join(chr_directory, chr_file)
        if chr_file.endswith('.gz'):
            gunzip_cmp = ('gunzip', chr_file)
            res0 = subprocess.Popen(gunzip_cmp)
            res0.wait()
            chr_file =  chr_file.rsplit('.', 1)[0]
            if res0.returncode != 0:
                exit(1)

        with open(chr_file) as f:
            for line in f:
                if line.startswith('>'):
                    chromo = line.strip().split()[0][1:]
                    if chromo in target_chromo_list:
                        chromo_file_dict[chromo] = chr_file
    return chromo_file_dict


def make_windows(df, temp_dir, step, output_file):
    # bedtools makewindows will lose information on strand
    # split df by strand
    # gtf is 1 based, bed is 0-based
    df['end'] += 1
    for strand in ['+', '-']:
        tmp_file = temp_dir + '.target.strand%s.bed' % strand
        df[df['strand'] == strand][['seqname', 'start', 'end']].sort_values(['seqname', 'start', 'end']).to_csv(
            tmp_file, sep='\t', header=False, index=False)
        bedtools_merge = ('bedtools', 'merge', '-i', tmp_file)
        bedtools_cmd = ('bedtools', 'makewindows', '-b', '-', '-w', '13', '-s', str(step))
        awk_cmd = ("awk", "{print $0" + '"\t"$1"_' + strand + '_"$2 "\t0\t' + strand + '"}')
        res0 = subprocess.Popen(bedtools_merge, stdout=subprocess.PIPE)
        res1 = subprocess.Popen(bedtools_cmd, stdin=res0.stdout, stdout=subprocess.PIPE)
        res2 = subprocess.Popen(awk_cmd, stdin=res1.stdout, stdout=open(output_file, 'a'))
        res2.wait()
        res1.wait()

        if res1.returncode == 0 and res2.returncode == 0:
            # clean up temp files
            os.remove(tmp_file)
        else:
            sys.exit(1)


def get_window_seq(chromo_dict, input_file, output_file):
    with open(input_file) as f:
        nb_cols = len(f.readline().split())
    if nb_cols == 0:
        os.rename(input_file, output_file)
    else:
        chromo_file = input_file + '.chromo'
        # gtf is 1-based and bed 0-based
        split_chromo = ("awk", ' $3 - $2 >= 13 {print $1 "\t" $2 - 1 "\t" $3 - 1 "\t"$'
                        + '"\t"$'.join(str(i) for i in range (4, nb_cols + 1))
                        + ' > ' + '"' + chromo_file + '"' + "$1" + "}", input_file)
        res0 = subprocess.Popen(split_chromo)
        res0.wait()
        if res0.returncode == 0:
            chromo_list = [f.replace(os.path.basename(input_file) + '.chromo', '') for f in os.listdir(os.path.dirname(chromo_file))
                           if os.path.basename(input_file) + '.chromo' in f]
            os.remove(input_file)
            for chromo in chromo_list:
                chr_fasta = chromo_dict[chromo]
                getfasta_cmd = ('bedtools', 'getfasta', '-fi', chr_fasta, '-bed', chromo_file + chromo, '-bedOut', '-s')
                awk_cmd = ('awk', '{print $1 "\t" $2 + 1 "\t" $3 + 1 "\t"$'
                        + '"\t"$'.join(str(i) for i in range (4, nb_cols + 2)) # +2 because we added sequence column
                        + '}' )
                res1 = subprocess.Popen(getfasta_cmd, stdout=subprocess.PIPE)
                res2 = subprocess.Popen(awk_cmd, stdin=res1.stdout, stdout=open(output_file, 'a'))
                res2.wait()
                res1.wait()
                if res1.returncode == 0 and res2.returncode == 0:
                    # clean up temp files
                    os.remove(chromo_file + chromo)
                else:
                    sys.exit(1)
        else:
            sys.exit(1)



def get_seq_bedtools(dataf, chr_directory, verbose, temp_dir, step):
    if verbose:
        print('Extracting target sequence')

    dataf['seqname'] = dataf['seqname'].map(str)
    target_chr_list = set(dataf['seqname'].unique())

    # list the chromosomes contained in each file
    chromo_file_dict = get_chr_file_dict(chr_directory, target_chr_list)

    missing_chromo = target_chr_list - set(chromo_file_dict.keys())
    if len(missing_chromo) > 0:
        print('Missing chromosome(s) in fasta files:\n%s\nExiting.' %(', '.join(missing_chromo)), file=sys.stderr)
        sys.exit(1)

    window_file = temp_dir + '.targetwindows.bed'
    window_seq_file = temp_dir + '.targetwindows.seq.bed'
    make_windows(dataf, temp_dir, step, window_file)
    get_window_seq(chromo_file_dict, window_file, window_seq_file)
    return chromo_file_dict
