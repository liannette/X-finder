import os
import subprocess
import sys
import zipfile


def check_for_stderr(output_stderr):
    """ 
    Checks if the subprocess had a stderr. If yes, the error message
    is printed and the program is stopped 
    """
    try:
        assert output_stderr == None
    except: 
        print(output_stderr)
        sys.exit(1)
    

def create_dir(outdir):
    """ Creates directory if it does not already exist"""
    if os.path.isdir(outdir) is False:
        os.mkdir(outdir)


def get_all_assemblies(genus, assembly_lvl):
    """
    Get all refseq assemblies of the genus, sorted by RefSeq assembly accession
    """
    args = [
        "ncbi-genome-download", "--dry-run", "--assembly-levels", assembly_lvl,
        "--genera", genus, "bacteria"
        ]
    output = subprocess.run(args, stdout=subprocess.PIPE, stdin=subprocess.PIPE, check=True)
        # check=True raises an CalledProcessError exception if process exits with a non-zero exit code
    check_for_stderr(output.stderr)
    assemblies = output.stdout.decode('utf-8').split("\n")[1:-1]
    return sorted(assemblies, reverse=True)


def filter_assemblies(assemblies):
    """
    Only one assembly for each strain, chooses the assembly with the highest accession.
    No assemblies with unknown species. 
    """

    filtered_assemblies = list()
    species_list = set(["sp."]) # exclude unknown species
    for assembly in assemblies:
        species = assembly.split()[2]
        if species not in species_list:
            filtered_assemblies.append(assembly)
            species_list.add(species)
    return filtered_assemblies


def refseq_accessions(genus, assembly_lvl, max_genomes, outdir):

    # Get all refseq assemblies of the genus
    assemblies = get_all_assemblies(genus, assembly_lvl)
    # only one assembly per organism name
    filtered_assemblies = filter_assemblies(assemblies)
    filtered_assemblies = filtered_assemblies[:max_genomes]
    # Write the assemblies to a file
    create_dir(outdir)
    filename = outdir / "assembly_accessions.txt"
    with open(filename, "w") as outfile:
        for assembly in filtered_assemblies:
            print(assembly.split()[0], file=outfile)
    return filename


def download_genomes(outdir, inputfile):
    
    # download
    outputfile = outdir / "ncbi_dataset.zip"
    args = [
        "datasets", "download", "genome", "accession", "--inputfile", 
        inputfile, "--include-gbff", "--exclude-gff3", 
        "--exclude-protein", "--exclude-rna", "--exclude-genomic-cds", 
        "--exclude-seq", "--filename", outputfile
        ]
    output = subprocess.run(args, stdout=subprocess.PIPE, stdin=subprocess.PIPE, check=True)
        # check=True raises CalledProcessError exception if process exits with a non-zero exit code
    check_for_stderr(output.stderr)
    sys.exit(1)
    # Unzip the file
    with zipfile.ZipFile(outputfile, 'r') as zip_ref:
        zip_ref.extractall(outdir)
    # delete unnecessary files
    os.remove(outputfile)
    os.remove(outdir / "README.md")  