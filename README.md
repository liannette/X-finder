# X-finder
X-finder is a command line tool for the detection of novel gene clusters 
for natural products through large scale genome comparison.

Developed by:

Annette Lien, Xinglin Jiang, Simon Shaw

E-mail: a.lien@posteo.de


## Dependencies

X-finder is build in python3.8. We recommend installing X-finder in a conda 
environment like so:

```
# navigate to where you want to download the xfinder repo
git clone https://github.com/liannette/X-finder.git

# navigate into the repo
cd X-finder

# create new environment
conda create -n xfinder python=3.8

# activate new environment
conda activate xfinder

# install X-finder and dependencies
pip install .
```

X-finder also requires diamond (version 2.0.14), which can be installed 
with conda in your xfinder environment:
```
conda install -c bioconda diamond=2.0.14
```

X-finder comes with an additional script to directly download 
genomes and run antiSMASH on them. This is the preferred way to generate
input genomes for X-finder. We recommend using antiSMASH version 6.1.1. More information on how to install antiSMASH: https://docs.antismash.secondarymetabolites.org/install/ 




## Usage

xfinder.py executes the main functionalities of X-finder. The required imput files for the X-finder analysis are two directories with genomes in GenBank format and one fasta file containing core genome coding sequences. A core genome file for Streptomyces is available in the data folder, but other core genomes can be generated with OrthoVenn2 (https://orthovenn2.bioinfotoolkits.net/home). The query and reference input genomes must contain annotations for coding sequences, PFAM domains, and antiSMASH cluster. We recommend to run antiSMASH (on minimal and fullhmmer settings) on each genome and then use the full genome .gbk file from the antiSMASH output.

The additional script get_genomes.py downloads genomes from NCBI and runs themthrough antiSMASH (on minimal and fullhmmer settings) in one step. get_genomes.py takes as input a file with one RefSeq accession per line.We suggest selecting genome assemblies (preferably at an assembly level of chromosome or complete) from https://www.ncbi.nlm.nih.gov/data-hub/genome/. A table of the selected assemblies can be downloaded and the RefSeq accessions extracted from the tsv file. 

### xfinder.py
```
usage: xfinder.py [-h] -o <dir> -q <dirs> -r <dirs> --core_genome <file> [--transporter_pfams <file>] [-t <int>] [--max_l50 <int>]
                  [--min_seed_size <int>] [--max_gap <int>] [--min_sublist_size <int>] [--min_dna_length <int>] [--core_cutoff <float>]
                  [--transporter_cutoff <float>]

optional arguments:
  -h, --help            show this help message and exit
  -o <dir>, --out_folder <dir>
                        Output directory, this will contain all output data files. Must not exist already.
  -q <dirs>, --queries_dirs <dirs>
                        Input directories of query gbk files, with antismash annotation (--fullhmmer). Multiple directories are possible,
                        paths must be seperated with a ','
  -r <dirs>, --ref_dirs <dirs>
                        Input directories of reference gbk files, with antismash annotation (--fullhmmer). Multiple directories are
                        possible, paths must be seperated with a ','
  --core_genome <file>  Fasta file containing core genome sequences.
  --transporter_pfams <file>
                        File containing all PFAM domains that are associated with transporter function, one PFAM number (without the PF)
                        per line. (default: None)
  -t <int>, --threads <int>
                        Set the number of threads the script may use (default: use all available cores)
  --max_l50 <int>       Only considers input genomes with a L50 value of this or lower. L50 is defined as count of smallest number of
                        contigs whose length sum makes up half of genome size (default: None)
  --min_seed_size <int>
                        The minimum number of pfam domains that match exactly during the initial step of finding common pfam domains
                        between a query and ref genome (each genome represented by a genomepfamsequence). Do not change this value,
                        except you know what you are doing! (default: 2)
  --max_gap <int>       The maximum number of pfam domains in a gap, when combining neighboring/overlapping exact matches (default: 2)
  --min_sublist_size <int>
                        Both sublists of a hit must contain at least this many PFAM domains. (default: 6)
  --min_dna_length <int>
                        The DNA sequence length of the query sublist must be at this long (default: 7000)
  --core_cutoff <float>
                        Cutoff for core genome filtering (default:0.5)
  --transporter_cutoff <float>
                        Cutoff for transporter function filtering (default:0.2)
```

Example command for xfinder.py:
```
python3 xfinder.py -o ../xfinder_output -q ../query_genomes -r ../ref_genomes1,../ref_genomes --core_genome data/Streptomyces_core_genome.fasta --transporter_pfams data/transporter_pfams.txt --max_l50 3 -t 8
```


### get_genomes.py

```
usage: get_genomes.py [-h] -o <dir> -i <str or file> [-t <int>] [-a <file>]

optional arguments:
  -h, --help            show this help message and exit
  -o <dir>, --output_dir <dir>
                        Output directory, this will contain all output data files.
  -i <str or file>, --refseq_accessions_path <str or file>
                        Absolute path to an input file containing one RefSeq accession per line.
  -t <int>, --threads <int>
                        Set the number of threads the script may use (default: use all available cores)
  -a <file>, --antismash_path <file>
                        Absolute path to the antismash program. Only necessary, if antismash is not loaded as a module and can't be
                        accessed by writing 'antismash' in the commandline.
```

Example command for get_genomes.py:
```
python3 get_genomes.py -o ../query_genomes -i ../query_genomes/refseq_acc.txt -t 8 -a /opt/antismash/6.1.1/antismash
```