"""
Author: Annette Lien (a.lien@posteo.de)

Additional script for simple generation of input files for X-finder

usage:
python xfinder.py -h
"""
import os
import subprocess
import sys
import glob
from subprocess import Popen, PIPE
import datetime
import argparse
from multiprocessing import cpu_count


def get_commands():
    parser = argparse.ArgumentParser(
        description= "Additional script to generate input files for X-finder. "
        "RefSeq assemblies are downloaded from NCBI and run through "
        "antiSMASH with --minimal and --fullhmmer setting. Takes as input a "
        "file with one RefSeq accession per line. We suggest getting the "
        "RefSeq accessions by downloading a tsv file from "
        "https://www.ncbi.nlm.nih.gov/data-hub/genome/ and then copying the "
        "RefSeq accessions into a text file.")
    parser.add_argument(
        "-o", "--output_dir", dest="outdir", required=True, metavar="<dir>",
        help="Output directory, this will contain all output data files.")
    parser.add_argument(
        "-i", "--refseq_accessions_path", dest="refseq_acc_file", 
        required=True, metavar="<file>", help="Path to an "
        "input file containing one RefSeq accession per line.")
    parser.add_argument(
        "-t", "--threads", dest="threads", default=cpu_count(),
        help="Set the number of threads the script may use (default: use all \
        available cores)", type=int, metavar="<int>")
    parser.add_argument(
        "-a", "--antismash_path", dest="antismash_path", default="antismash",
        metavar="<file>", help="Path to antismash. Only necessary, if "
        "antismash is not loaded as a module and can't be accessed by writing "
        "'antismash' in the commandline.")
    return parser.parse_args()


def print_log(msg, out_dir=None):
    ''' Adds a timestamp to a string and prints it to log and stdout'''
    now = datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    if out_dir is not None:
        with open(os.path.join(out_dir, "log.txt"), "a") as outfile:
            print(f"{now}\t{msg}", file=outfile)
    print(f"{now}\t{msg}")
    sys.stdout.flush()


def _print_stdout_stderr(stdout, stderr, outdir):
    """ Prints stdout and stderr of a subprocess to log and stdout"""
    for msg in (stdout, stderr):
        if len(msg) > 0:
            print_log(msg.decode(), outdir)


def download_genomes(refseq_acc_file, outdir):

    if os.path.isdir(outdir) is False:
        os.mkdir(outdir)
    if os.path.exists(refseq_acc_file) is False:
        raise RuntimeError("The assembly accession file does not exist: "
                           f"{refseq_acc_file}")
    
    # Dry run to count the number of genomes to download
    args = ["ncbi-genome-download", "--dry-run", "--assembly-accessions", 
            refseq_acc_file, "bacteria"]
    process = Popen(args, stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate() 
    if process.returncode != 0:   
        _print_stdout_stderr(stdout, stderr, outdir)
        raise RuntimeError("Error while running ncbi-genome-download")     
    num_assemblies = len(stdout.decode().split("\n")[1:-1])
    print_log(f"Downloading {num_assemblies} assemblies", outdir)

    # Actually downloading the genomes
    args = [
        "ncbi-genome-download", "--human-readable", 
        "--assembly-accessions", refseq_acc_file,
        "--output-folder", outdir,
        "bacteria"
        ]
    process = Popen(args, stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate() 
    _print_stdout_stderr(stdout, stderr, outdir)


def _get_gbff_files(outdir):
    # Get all genomes files
    genomic_files = glob.glob(
        os.path.join(outdir, "refseq", "bacteria", '*', '*.gbff.gz')
        )
    return genomic_files

def _run_antismash(gbff_file, outdir, outdir_antismash, threads, antismash_path):
    output_basename = os.path.basename(gbff_file).rstrip(".gbff.gz")
    args = [antismash_path, "--minimal", "--fullhmmer", "--output-dir", 
            os.path.join(outdir_antismash, output_basename),"--cpus", 
            str(threads), gbff_file,]
    process = Popen(args, stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate() 
    _print_stdout_stderr(stdout, stderr, outdir)


def _create_links(outdir):
    """ Create the symbolic links for the annotated genomes files """
    
    # Destination directory for all annotated files
    dst_dir = os.path.join(outdir, "final")
    if os.path.isdir(dst_dir) is False:
        os.mkdir(dst_dir)
    
    for dir_name in glob.glob(os.path.join(outdir, "antismash", "*/")):
        genome_gbk = dir_name + dir_name.split("/")[-2] + ".gbk"
        if os.path.isfile(genome_gbk):
            # Create the link
            dst_path = os.path.join(dst_dir, os.path.basename(genome_gbk))
            process = subprocess.Popen(["ln", "-s", genome_gbk, dst_path], 
                                       stdout=PIPE, stderr=PIPE)
            stdout, stderr = process.communicate() 
            _print_stdout_stderr(stdout, stderr)


def run_antismash(outdir, threads, antismash_path):
    """ Runs antiSMASH on all file downloaded previouly with 
    ncbi-genome-download. Settings for antiSMASH are --minimal and 
    --fullhmmer. Creates a directory with the name 'final' in the output
    directory, containing symbolic links to the genome genbank files 
    from the antiSMASH output."""

    # Get a list of all paths of the genome files
    genomic_files = _get_gbff_files(outdir)
    
    # Create the destination directory
    outdir_antismash = os.path.join(outdir, "antismash")
    if os.path.isdir(outdir_antismash) is False:
        os.mkdir(outdir_antismash)
    print_log(f"Writing antismash results to dir: {outdir_antismash}", outdir)
    
    # Run antismash on commandline
    for i in range(len(genomic_files)):
        gbff_file = genomic_files[i]
        print_log(f"Start antismash on file {i+1}/{len(genomic_files)}: "
                  f"{os.path.basename(gbff_file)}", outdir)
        _run_antismash(gbff_file, outdir, outdir_antismash, threads, 
                   antismash_path)
        
    # symbolic links for the antismash genbank files   
    _create_links(outdir)
    
    
if __name__ == "__main__":
    
    args = get_commands()
    
    download_genomes(args.refseq_acc_file, args.outdir)
    run_antismash(args.outdir, args.threads, args.antismash_path)