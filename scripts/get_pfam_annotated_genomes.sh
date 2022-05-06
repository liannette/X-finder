#!/bin/bash

GENUS=bacillus
ASSEMBLY_LEVEL=complete,chromosome # Assembly levels of genomes to download. A comma-separated list of assembly levels is also possible. For example: "complete,chromosome". Choose from: ['all', 'complete', 'chromosome', 'scaffold', 'contig']
MAX_GENOMES=100
CPUS=8 # number of CPUs for antismash

SCRIPT_DIR=$PWD
OUTDIR="$( dirname "$PWD" )/data/genomes/$GENUS"
mkdir $OUTDIR
> $OUTDIR/output.log



# Download genomes, only one strain of each organism

python3 $SCRIPT_DIR/download_genomes.py --genus $GENUS --assembly-levels $ASSEMBLY_LEVEL --max_genomes $MAX_GENOMES --outdir $OUTDIR \
    >> $OUTDIR/output.log 2>&1

unzip $OUTDIR/ncbi_dataset.zip -d $OUTDIR
rm $OUTDIR/ncbi_dataset.zip
rm $OUTDIR/README.md



# PFAM annotate the genomes with antismash

python3 $SCRIPT_DIR/create_symbolic_links.py --mode 1 --outdir $OUTDIR \
    >> $OUTDIR/output.log 2>&1

cd $OUTDIR/antismash
module load antismash
ls | xargs -L 1 -d '\n' antismash --minimal --fullhmmer --cpus $CPUS \
    >> $OUTDIR/output.log 2>&1
rm *.gbff


# Create the symbolic links for discoverBGC

python3 $SCRIPT_DIR/create_symbolic_links.py --mode 2 --outdir $OUTDIR \
    >> $OUTDIR/output.log 2>&1
