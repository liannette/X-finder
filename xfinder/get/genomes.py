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


def get_all_assemblies(genus, tax_ids, assembly_acc, assembly_lvl):
    """
    Get all refseq assemblies of the genus, sorted by RefSeq assembly accession
    """
    
    assert (genus, tax_ids, assembly_acc).count(None) == 2
    
    args = ["ncbi-genome-download", "--dry-run", 
            "--assembly-levels", assembly_lvl]
    if genus is not None:
        args += ["--genera", genus]
    elif tax_ids is not None:
        args += ["--taxids", tax_ids]
    else:
        args += ["--assembly-accessions", assembly_acc]
    args.append("bacteria")

    # check=True raises an CalledProcessError exception if process 
    # exits with a non-zero exit code
    output = subprocess.run(args, stdout=subprocess.PIPE, 
                            stdin=subprocess.PIPE, check=True)
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


def get_refseq_accessions(genus, tax_id, assembly_acc, assembly_lvl, 
                          max_genomes, one_per_species, outdir):

    # Get all refseq assemblies of the genus
    assemblies = get_all_assemblies(genus, tax_id, assembly_acc, assembly_lvl)
    # only one assembly per organism name
    if one_per_species:
        assemblies = filter_assemblies(assemblies)
    if max_genomes is not None:
        assemblies = assemblies[:max_genomes]
        
    # Write the assemblies to a file
    if os.path.isdir(outdir) is False:
        os.mkdir(outdir)
    filename = os.path.join(outdir, "refseq_assembly_accessions.txt")
    with open(filename, "w") as outfile:
        for assembly in assemblies:
            print(assembly.split()[0], file=outfile)
    return filename


def download_genomes(outdir, inputfile):
    
    outputfile = os.path.join(outdir, "ncbi_dataset.zip")
    args = [
        "datasets", "download", "genome", "accession", "--inputfile", 
        inputfile, "--include-gbff", "--exclude-gff3", 
        "--exclude-protein", "--exclude-rna", "--exclude-genomic-cds", 
        "--exclude-seq", "--filename", outputfile
        ]
    # check=True raises CalledProcessError exception if process exits 
    # with a non-zero exit code
    output = subprocess.run(args, stdout=subprocess.PIPE, 
                            stdin=subprocess.PIPE, check=True)
    check_for_stderr(output.stderr)
    # Unzip the file
    with zipfile.ZipFile(outputfile, 'r') as zip_ref:
        zip_ref.extractall(outdir)
    # delete unnecessary files
    os.remove(outputfile)
    os.remove(os.path.join(outdir,"README.md"))