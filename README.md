# **snoGloBe**: a box C/D snoRNA interaction predictor

## **Getting Started**

snoGloBe is supported on Linux (tested on Ubuntu, Fedora, CentOS).

### **Dependencies**

snoGloBe works with python 3.5 or greater and requires the following packages available from pip:

* pandas
* biopython
* scikit-learn v0.21.3

snoGloBe also needs 

* [bedtools](http://bedtools.readthedocs.io/en/latest/). 
[Click here for installation guidelines](http://bedtools.readthedocs.io/en/latest/content/installation.html)

### **Installing**


First, you must clone the git repository:
```
cd /path/to/clone/
git clone http://gitlabscottgroup.med.usherbrooke.ca/scott-group/snoglobe.git
```

For easier use, add the /bin/ path of the snoglobe repository to your PATH environment variable.
This is easily done by opening your .bashrc file in your $HOME directory and adding this line at the bottom:

```
export PATH=/path/to/clone/snoglobe/bin:$PATH
```

You should be able to execute the **snoglobe** command from your terminal.

 ```
snoglobe --help
usage: snoglobe.py [-h] [-v] [-n NB_THREADS] [-s STEPSIZE] [-c CHUNKSIZE]
                   [-t THRESHOLD] [-m] [-w NB_WINDOWS] [--seq] [--verbose]
                   sno_fasta target_ids gtf chromo_fasta_dir output

positional arguments:
  sno_fasta             fasta file containing snoRNA sequences
  target_ids            txt file containing target identifiers (gene_id,
                        transcript_id or exon_id), ids should match the ones
                        in the gtf file
  gtf                   Annotation file in .gtf format. Preferably an
                        annotation of the whole genome or whole chromosomes of
                        specified targets
  chromo_fasta_dir      Directory containing fasta files of individual
                        chromosome
  output                Output file name

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         Show the version and exit
  -n NB_THREADS, --nb_threads NB_THREADS
                        Number of threads to be used. Default: 1
  -s STEPSIZE, --stepsize STEPSIZE
                        Number of nucleotides between each target sliding
                        window. Default: 1
  -c CHUNKSIZE, --chunksize CHUNKSIZE
                        Number of combinations to run at once. Default:
                        2500000
  -t THRESHOLD, --threshold THRESHOLD
                        Minimum score for an interaction to be reported, value
                        between 0 and 1. Default: 0.95
  -m, --merge           Use this option to merge consecutive positive windows
                        as one
  -w NB_WINDOWS, --nb_windows NB_WINDOWS
                        Minimum number of consecutive windows for an
                        interaction to be reported. Must be used with -m
                        option. Default: 1
  --seq                 Add target and snoRNA interaction sequences to the
                        output
  --verbose             Print the steps
```

## **Recommendation**

Running a single snoRNA against every human genes will take a few days. We recommend focusing on a subset of genes of
interest and using a cluster for more exhaustive jobs.

## **Authors**

* **Gabrielle Deschamps-Francoeur**
* **Michelle S Scott** - [website](https://bioinfo-scottgroup.med.usherbrooke.ca/)

Questions should be directed to: _gabrielle.deschamps-francoeur@usherbrooke.ca_


## **License**

snoGloBe : a box C/D snoRNA interaction predictor
Copyright (C) 2020 Gabrielle Deschamps-Francoeur & Michelle S Scott

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

For more information about the GNU General Public License, see <http://www.gnu.org/licenses/>.


## **How to cite**

Please cite as :

Deschamps-Francoeur G, Couture S, Abou-Elela S, Scott MS. The snoGloBe interaction predictor enables a broader study 
of box C/D snoRNA functions and mechanisms. bioRxiv. 2021. Available from: 
<https://www.biorxiv.org/content/10.1101/2021.09.14.460265v1>