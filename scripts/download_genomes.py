import os
import subprocess
import sys
import argparse
from datetime import datetime


def dry_run_ncbi_genome_download(genus, assembly_lvl):
    args = [
        "ncbi-genome-download", "--dry-run", "--assembly-levels", assembly_lvl,
        "--genera", genus, "bacteria"
    ]
    output = subprocess.run(args, stdout=subprocess.PIPE, stdin=subprocess.PIPE, check=True)
    try:
        assert output.stderr == None
    except: 
        print(output.stderr)
        sys.exit(1)
    assemblies = output.stdout.decode('utf-8').split("\n")[1:-1]

    return output.stdout.decode('utf-8')


def filter_assemblies(assemblies, max_assemblies):
    """
    Outputs only one assembly per organism name. Chooses the first assembly in the list.
    """

    filtered_assemblies = list()
    assembly_accessions = list()
    organisms = set()
    for assembly in assemblies:
        assembly_accession, organism_name, strain = assembly.split("\t")
        if organism_name not in organisms:
            filtered_assemblies.append(assembly)
            assembly_accessions.append(assembly_accession)
            organisms.add(organism_name)

    return filtered_assemblies[:max_assemblies], assembly_accessions[:max_assemblies]
        

def dir_name():
    DIR = os.path.dirname(os.path.abspath(__file__))
    return DIR


def download_genomes(genus, accessions, outdir):
    args = [
        "datasets", "download", "genome", "accession", *accessions, 
        "--include-gbff", "--exclude-gff3", "--exclude-protein", 
        "--exclude-rna", "--exclude-genomic-cds", "--exclude-seq",
        "--filename", "{}/ncbi_dataset.zip".format(outdir)
    ]
    output = subprocess.run(args, stdout=subprocess.PIPE, stdin=subprocess.PIPE, check=True)
    try:
        assert output.stderr == None
    except: 
        print(output.stderr)
        sys.exit(1)
    
    

def main(genus, assembly_lvl, max_genomes, outdir):

    # Get all refseq assemblies of the genus
    stdout = dry_run_ncbi_genome_download(genus, assembly_lvl)
    assemblies = stdout.split("\n")[1:-1]

    # only one assembly per organism name
    assemblies, assembly_accessions = filter_assemblies(assemblies, max_assemblies=max_genomes)

    # Write the assembly accession and organism names to file
    with open(os.path.join(outdir, genus+"_download.txt"), "w") as infile:
        infile.write("\n".join(assemblies))
        
    # Download the files
    download_genomes(genus, assembly_accessions, outdir)


if __name__ == "__main__":

    print("Commandline input:", " ".join(sys.argv))

    parser = argparse.ArgumentParser()
    parser.add_argument('--genus', action="store", dest="genus", type=str, help='download genomes of this genus name')
    parser.add_argument('--assembly-levels', action="store", dest="assembly_lvl", type=str, help='Assembly levels of genomes to download. A comma-separated list of assembly levels is also possible. For example: "complete,chromosome". Choose from: ["all", "complete", "chromosome", "scaffold", "contig"]')
    parser.add_argument('--max_genomes', action="store", dest="max_genomes", type=int, help='Will not download more genomes than this number')
    parser.add_argument('--outdir', action="store", dest="outdir", type=str, help='Directory to save the genomes')
    
    args = parser.parse_args()
    genus = args.genus
    assembly_lvl = args.assembly_lvl
    max_genomes = args.max_genomes
    outdir = args.outdir

    print("{}: Program started.".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
    main(genus, assembly_lvl, max_genomes, outdir)
    print("{}: Program finished".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))