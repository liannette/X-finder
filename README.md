# X-finder
X-finder is a command line tool for the detection of novel gene clusters for natural products through large scale genome comparison.

Developed by:

Xinglin Jiang
Annette Lien
Simon Shaw

E-mail: a.lien@posteo.de


## Dependencies

X-finder is build and tested in python3.8. We recommend installing
X-finder in a conda environment like so:

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

X-finder also requires diamond (version 2.0.14), which can be installed with conda in your xfinder environment:
```
conda install -c bioconda diamond=2.0.14
```


## Usage

All paths must be absolute.

```
usage: xfinder.py [-h] -o <dir> -q <dir>[,<dir>] -r <dir>[,<dir>] --core_genome <file> [--transporter_pfams <file>] [-t <int>] [--max_l50 <int>]
                  [--min_seed_size <int>] [--max_gap <int>] [--min_sublist_size <int>] [--min_dna_length <int>] [--core_cutoff <float>]
                  [--transporter_cutoff <float>]

optional arguments:
  -h, --help            show this help message and exit
  -o <dir>, --out_folder <dir>
                        Output directory, this will contain all output data files. Must not exist already.
  -q <dir>[,<dir>], --queries_dirs <dir>[,<dir>]
                        Input directories of query gbk files. Multiple directories possible, seperate paths with a ','
  -r <dir>[,<dir>], --ref_dirs <dir>[,<dir>]
                        Input directories of query gbk files. Multiple directories possible, seperate paths with a ','
  --core_genome <file>  Fasta file containing core genome sequences.
  --transporter_pfams <file>
                        File containing all PFAM domains that are associated with transporter function, one PFAM number (without the PF) per line.
                        (default: None)
  -t <int>, --threads <int>
                        Set the number of threads the script may use (default: use all available cores)
  --max_l50 <int>       Only considers input genomes with a L50 value of this or lower. L50 is defined as count of smallest number of contigs whose
                        length sum makes up half of genome size (default: None)
  --min_seed_size <int>
                        The minimum number of pfam domains that match exactly during the initial step of finding common pfam domains between a query and
                        ref genome (each genome represented by a genomepfamsequence). Do not change this value, except you know what you are doing!
                        (default: 2)
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

Example command:
```
python3 xfinder.py -o OUTDIR -q  -r /data/s202633/X-finder_inputs/genomes/ref_genomes --core_genome /data/s202633/X-finder_inputs/core_genome/Streptomyces.fasta --transporter_pfams /data/s202633/X-finder_inputs/transporter_pfams/all_transporter_pfams.txt --max_l50 3 -t 8
```