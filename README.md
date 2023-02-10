# X-finder
X-finder is a command line tool for the detection of novel gene clusters 
for natural products through large scale genome comparison.

Developed by:

Annette Lien, Xinglin Jiang, Simon Shaw

E-mail: a.lien@posteo.de


## Dependencies

X-finder is build in python3.8. We recommend installing X-finder in a conda environment like so:

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


If there are problems using the above descriped installation method, one can try the following instead:
```
# navigate to where you want to download the xfinder repo
git clone https://github.com/liannette/X-finder.git

# navigate into the repo
cd X-finder

# create new conda environment and install requirements in one go
conda env create -f conda_env.yml

# activate new environment
conda activate xfinder

# install X-finder
pip install .
```

X-finder comes with an additional script to directly download 
genomes and run antiSMASH on them. This is the preferred way to generate input genomes for X-finder. We recommend using antiSMASH version 6.1.1. More information on how to install antiSMASH: https://docs.antismash.secondarymetabolites.org/install/ 




## Usage

X-finder contains two commandline scripts:
* xfinder.py for main functionalities
* getgenomes.py for generating input genomes for xfinder


### xfinder.py

xfinder.py executes the main functionalities of X-finder. The required imput files are at least two directories with genomes in GenBank format (one with reference genomes and one with query genomes) and one fasta file containing coding sequences (CDS) of the core genomes for filtering the results. 

The query and reference input genomes must contain annotations for CDS and PFAM domains. Ideally they also contain antiSMASH cluster annotation. We recommend to run antiSMASH (on minimal and fullhmmer settings) on each genome and then use the full genome .gbk file from the antiSMASH output. This can be automated by using getgenomes.py.

xfinder.py also needs a core genome file as input. A Streptomyces core genome file is available in the data folder, but other core genomes can be generated with OrthoVenn2 (https://orthovenn2.bioinfotoolkits.net/home). 

We also provide a file containing the transporter PFAM numbers in the data folder. However, if no file is specified, xfinder.py will generate its own list of transporter PFAM numbers. This prolongs the runtime by approximately 15 minutes.

Toggle -h or --help for additional command line arguments and default values: 
```
# Show help message
python3 xfinder.py -h
```
Help message:
```
usage: xfinder.py [-h] -o <dir> -q <dirs> -r <dirs> --core_genome <file>
                  [--transporter_pfams <file>] [-t <int>] [--max_l50 <int>]
                  [--min_seed_size <int>] [--max_gap <int>] [--min_sublist_size <int>]
                  [--min_dna_length <int>] [--core_cutoff <float>]
                  [--transporter_cutoff <float>]

optional arguments:
  -h, --help            show this help message and exit
  -o <dir>, --output_dir <dir>
                        Output directory, this will contain all output data files. Must
                        not exist already.
  -q <dirs>, --queries_dirs <dirs>
                        Input directories of query gbk files, Genomes must be in
                        Genbank format with '.gbk' file extension and contain antismash
                        annotation (--fullhmmer). Multiple directories are possible,
                        paths must be seperated with a ','
  -r <dirs>, --ref_dirs <dirs>
                        Input directories of reference genomes. Genomes must be in
                        Genbank format with '.gbk' file extension and contain antismash
                        annotation (--fullhmmer). Multiple directories are possible,
                        paths must be seperated with a ','
  --core_genome <file>  Fasta file containing core genome sequences.
  --transporter_pfams <file>
                        File containing all PFAM domains that are associated with
                        transporter function, one PFAM number (without the PF) per
                        line. (default: None)
  -t <int>, --threads <int>
                        Set the number of threads the script may use (default: use all
                        available cores)
  --max_l50 <int>       Only considers input genomes with a L50 value of this or lower.
                        L50 is defined as count of smallest number of contigs whose
                        length sum makes up half of genome size (default: None)
  --min_seed_size <int>
                        The minimum number of pfam domains that match exactly during
                        the initial step of finding common pfam domains between a query
                        and ref genome (each genome represented by a
                        genomepfamsequence). Do not change this value, except you know
                        what you are doing! (default: 2)
  --max_gap <int>       The maximum number of pfam domains in a gap, when combining
                        neighboring/overlapping exact matches (default: 2)
  --min_sublist_size <int>
                        Both sublists of a hit must contain at least this many PFAM
                        domains. (default: 6)
  --min_dna_length <int>
                        The DNA sequence length of the query sublist must be at this
                        long (default: 7000)
  --core_cutoff <float>
                        Cutoff for core genome filtering (default: 0.5)
  --transporter_cutoff <float>
                        Cutoff for transporter function filtering (default: 0.2)
```

Example commands for xfinder.py:
```
# Only one directory for reference genomes and query genomes
python3 xfinder.py \
    --output_dir ../xfinder_output \
    --queries_dirs ../query_genomes \
    --ref_dirs ../ref_genomes \
    --max_l50 3 \
    --core_genome data/Streptomyces_core_genome.fasta \
    --transporter_pfams data/transporter_pfams.txt \
    --threads 8

# Several directories for reference genomes and query genomes
python3 xfinder.py \
    --output_dir ../xfinder_output \
    --queries_dirs ../query_genomes1,../query_genomes2 \
    --ref_dirs ../ref_genomes1,../ref_genomes2 \
    --max_l50 3 \
    --core_genome data/Streptomyces_core_genome.fasta \
    --transporter_pfams data/transporter_pfams.txt \
    --threads 8
```


### getgenomes.py

The additional script getgenomes.py downloads genomes from NCBI and runs them through antiSMASH (on minimal and fullhmmer settings). The resulting files can be used as input genomes for xfinder.py and are located unter <output_dir>/antismash/all. getgenomes.py takes as input a file containing RefSeq accessions, one per line. To generate the input file, we suggest selecting genome assemblies (preferably at an assembly level 'chromosome' or 'complete') from https://www.ncbi.nlm.nih.gov/data-hub/genome/. A table of the selected assemblies can be downloaded from that site and the RefSeq accessions extracted from that tsv file. 

Downloading the genomes and running antiSMASH can be done in two seperate steps or in one. You can choose between 'download' (only download), 'antismash' (only antiSMASH) or 'complete' (download and antiSMASH). Toggle -h or --help for additional command line arguments and default values: 

```
# Show help messages
python3 getgenomes.py -h
python3 getgenomes.py download -h
python3 getgenomes.py antismash -h
python3 getgenomes.py complete -h
```

Example commands for getgenomes.py:
```
# Download refseq assemblies from ncbi
python3 getgenomes.py download \
    --output_dir ../genomes \
    --refseq_accessions_path ../genomes/refseq_acc.txt \

# Run antiSMASH on the previously downloaded genomes
python3 getgenomes.py antismash \
    --output_dir ../genomes \
    --input_dir ../genomes/ncbi-download/all \
    --threads 8 \
    --antismash_path /opt/antismash/6.1.1/antismash

# Download genomes and run antiSMASH in one step
python3 getgenomes.py complete \
    --output_dir ../genomes \
    --refseq_accessions_path ../genomes/refseq_acc.txt \
    --threads 8 \
    --antismash_path /opt/antismash/6.1.1/antismash
```