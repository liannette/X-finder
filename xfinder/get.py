from pathlib import Path
import get.genomes
import get.annotation

# Variables
genus = "pseudomonas"
assembly_lvl = "complete"
max_genomes = 10
outdir = Path(__file__).resolve().parent.parent / "data" / "genomes" / genus
cpus=8

# Get a list of refseq assembly accessions for download
acc_file = get.genomes.refseq_accessions(genus, assembly_lvl, max_genomes, outdir)

# Download assemblies from ncbi
get.genomes.download_genomes(outdir, acc_file)
# Add PFAM Annotation using antismash
#get.annotation.main(outdir, cpus)