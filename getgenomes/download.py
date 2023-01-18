from subprocess import Popen, PIPE
import os
from getgenomes.common import print_log, print_stdout_stderr, create_dir
import glob


def _get_genome_files(outdir):
    # Get all genomes files
    genome_files = glob.glob(
        os.path.join(outdir, "refseq", "bacteria", '*', '*.gbff.gz')
        )
    return genome_files


def _create_links(outdir_downloads, outdir):
    """ Create the symbolic links for the downloaded genomes files """
    
    # Destination directory for all downloaded files
    dst_dir = os.path.join(outdir_downloads, "all")
    create_dir(dst_dir)
    for genome_file in _get_genome_files(outdir_downloads):
        process = Popen(["ln", genome_file, dst_dir], stdout=PIPE, stderr=PIPE)
        stdout, stderr = process.communicate() 
        print_stdout_stderr(stdout, stderr, outdir)     
    return dst_dir


def _count_genomes(refseq_acc_file, outdir):
    # Dry run to count the number of genomes to download
    args = ["ncbi-genome-download", 
            "--dry-run", 
            "--assembly-accessions", str(refseq_acc_file), 
            "bacteria"]
    process = Popen(args, stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate() 
    if process.returncode != 0:   
        print_stdout_stderr(stdout, stderr, outdir)
        msg = ("Error occured while running the following command: "
               f"{' '.join(args)}")
        print_log(msg, outdir)
        raise RuntimeError(msg)     
    num_assemblies = len(stdout.decode().split("\n")[1:-1])
    print_log(f"Downloading {num_assemblies} genomes", outdir)


def download_genomes(refseq_acc_file, outdir):
    
    # Check the arguments
    if os.path.isfile(refseq_acc_file) is False:
        raise RuntimeError(f"Not a file: {refseq_acc_file}")
    create_dir(outdir)
    
    # Dry run to count the number of genomes to download
    _count_genomes(refseq_acc_file, outdir)
    
    # download the genomes
    outdir_downloads = os.path.join(outdir, "ncbi-download")
    create_dir(outdir_downloads)
    args = [
        "ncbi-genome-download", 
        "--assembly-accessions", str(refseq_acc_file),
        "--output-folder", str(outdir_downloads),
        "bacteria"
        ]
    process = Popen(args, stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate() 
    print_stdout_stderr(stdout, stderr, outdir)
    if process.returncode != 0:
        msg = ("Error occured while running the following command: "
               f"{' '.join(args)}")
        print_log(msg, outdir)
        raise RuntimeError(msg)   
    
    # Create a directory with symbolic links of all downloaded genomes
    final_dir = _create_links(outdir_downloads, outdir)
    
    return final_dir