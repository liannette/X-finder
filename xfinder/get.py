#from pathlib import Path
import get.genomes
import get.annotation

### Choose one of the below

# Only download sequences of the provided genera. A comma-seperated 
# list of genera is also possible. For example: 
# "Streptomyces coelicolor,Escherichia coli"
genus =  None # "myxococcales"
# Only download sequences of the provided NCBI taxonomy IDs. A 
# comma-separated list of taxids is also possible. For example: 
# "9606,9685"
tax_id = None
# Only download sequences matching the provided NCBI assembly 
# accession(s). A comma-separated list of accessions is possible, as 
# well as a path to a filename containing one accession per line.
assembly_acc = "/data/s202633/X-finder_inputs/genomes/myxococcales_all/refseq_assembly_accessions.txt"

### More options

# Assembly levels of genomes to download. A comma-separated list of 
# assembly levels is also possible. For example: "complete,chromosome". 
# Choose from: ['all', 'complete', 'chromosome', 'scaffold', 'contig']
assembly_lvl = "all"
# Limit on to how many genomes to download. Will download the genomes 
# with the highest refseq accession
max_number_of_genomes = None
# Whether to download only one genome of each species. Will download the  
# genome with the highest refseq accession
one_per_species = False 
# Where to download the genomes to
outdir = "/data/s202633/X-finder_inputs/genomes/myxococcales_all"  #Path(__file__).resolve().parent.parent / "data" / "genomes" / genus
# How many cpus antismash will use
cpus = 8


# Get a list of refseq assembly accessions for download
# acc_file = get.genomes.get_refseq_accessions(
#     genus=genus, 
#     tax_id=tax_id, 
#     assembly_acc=assembly_acc,
#     assembly_lvl=assembly_lvl,
#     max_genomes=max_number_of_genomes, 
#     one_per_species=one_per_species,
#     outdir=outdir,
#     )

# # Download assemblies from ncbi
# get.genomes.download_genomes(outdir, acc_file)
# Add PFAM Annotation using antismash
get.annotation.main(outdir, cpus)