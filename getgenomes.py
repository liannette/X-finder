"""
Author: Annette Lien (a.lien@posteo.de)

Additional script for easy generation of input files for X-finder

usage:
python3 getgenomes.py -h
python3 getgenomes.py complete -h
python3 getgenomes.py download -h
python3 getgenomes.py antismash -h
"""

import argparse
import os
import sys
from getgenomes.common import create_dir, print_log
from getgenomes.download import download_genomes
from getgenomes.antismash import run_antismash


def print_help(commands):
    parser = argparse.ArgumentParser(
            description= "Script to generate input genomes for xfinder.py "
            "Please specify the mode: 'download' to download RefSeq genomes, "
            "'antismash' to do a minimal run of antiSMASH including full "
            "whole-genome HMMer analysis and 'complete' to do both in one "
            "step.", prog='python3 ' + os.path.basename(__file__))
    parser.add_argument('command', metavar="<command>", 
                        help=f"Choose from: {', '.join(commands)}")
    parser.print_help()


def get_commands_complete(prog_name):
    parser = argparse.ArgumentParser(
        description= "Script to generate input genomes for xfinder.py "
        "Takes as input a file with one RefSeq accession per line. Downloads "
        "the RefSeq assemblies from NCBI and runs them through antiSMASH with "
        "--minimal and --fullhmmer setting. We suggest getting the RefSeq "
        "accessions by downloading a tsv file from "
        "https://www.ncbi.nlm.nih.gov/data-hub/genome/ and then copying the "
        "RefSeq accessions into a text file.", prog=prog_name)
    parser.add_argument(
        "-o", "--output_dir", dest="outdir", required=True, metavar="<dir>",
        help="Output directory, this will contain all output data files.")
    parser.add_argument(
        "-i", "--refseq_accessions_path", dest="refseq_acc_file", 
        required=True, metavar="<file>", help="Path to an "
        "input file containing one RefSeq accession per line.")
    parser.add_argument(
        "-t", "--threads", dest="threads", default=16,
        help="Set the number of threads the antiSMASH may use, max value is "
        "16 (default: 16)", type=int, metavar="<int>")
    parser.add_argument(
        "-a", "--antismash_path", dest="antismash_path", default="antismash",
        metavar="<file>", help="Path to antismash. Only necessary, if "
        "antismash is not loaded as a module and can't be accessed by writing "
        "'antismash' in the commandline.")
    parser.add_argument('mode', help=argparse.SUPPRESS)
    return parser.parse_args()


def get_commands_download(prog_name):
    parser = argparse.ArgumentParser(
        description= "Script to generate input genomes for xfinder.py "
        "Takes as input a file with one RefSeq accession per line. The "
        "RefSeq assemblies are downloaded from NCBI. We suggest getting "
        "the RefSeq accessions by downloading a tsv file from "
        "https://www.ncbi.nlm.nih.gov/data-hub/genome/ and then copying the "
        "RefSeq accessions into a text file.", prog=prog_name)
    parser.add_argument('mode', help=argparse.SUPPRESS)
    parser.add_argument(
        "-o", "--output_dir", dest="outdir", required=True, metavar="<dir>",
        help="Output directory, this will contain all output data files.")
    parser.add_argument(
        "-i", "--refseq_accessions_path", dest="refseq_acc_file", 
        required=True, metavar="<file>", help="Path to an "
        "input file containing one RefSeq accession per line.")
    return parser.parse_args()


def get_commands_antismash(prog_name):
    parser = argparse.ArgumentParser( 
        description= "Script to generate input genomes for xfinder.py "
        "Run antiSMASH on all files of the input directory. The ideal "
        "input for antiSMASH is an annotated nucleotide file in Genbank "
        "format or EMBL format. antiSMASH is run with --minimal and "
        "--fullhmmer settings", prog=prog_name)
    parser.add_argument('mode', help=argparse.SUPPRESS)
    parser.add_argument(
        "-o", "--output_dir", dest="outdir", required=True, metavar="<dir>",
        help="Output directory, this will contain all output data files.")
    parser.add_argument(
        "-i", "--input_dir", dest="indir", required=True, metavar="<dir>", 
        help="Path to a directory containing genome files which will be used "
        "as input for antiSMASH.")
    parser.add_argument(
        "-t", "--threads", dest="threads", default=16,
        help="Set the number of threads the antiSMASH may use, max value is "
        "16 (default: 16)", type=int, metavar="<int>")
    parser.add_argument(
        "-a", "--antismash_path", dest="antismash_path", default="antismash",
        metavar="<file>", help="Path to antismash. Only necessary, if "
        "antismash is not loaded as a module and can't be accessed by writing "
        "'antismash' in the commandline.")
    return parser.parse_args()


def _prog_name(mode):
    return 'python3 ' + os.path.basename(__file__) + " " + mode


def command2log(outdir):
    create_dir(outdir)
    print_log(f"Command: {' '.join(sys.argv)}")


def download_and_antismash(mode):
    args = get_commands_complete(_prog_name(mode))
    command2log(args.outdir)
    genomes_dir = download_genomes(
        args.refseq_acc_file, 
        args.outdir,
        )
    run_antismash(
        genomes_dir, 
        args.outdir, 
        args.threads,
        args.antismash_path,
        )


def only_download(mode):
    args = get_commands_download(_prog_name(mode))
    command2log(args.outdir)
    download_genomes(
        args.refseq_acc_file, 
        args.outdir,
        )
    
    
def only_antismash(mode):
    args = get_commands_antismash(_prog_name(mode))
    command2log(args.outdir)
    run_antismash(
        args.indir, 
        args.outdir, 
        args.threads,
        args.antismash_path,
        )


if __name__ == '__main__':
    
    commands = ["complete", "download", "antismash"]

    # only the program name without additional arguments
    if len(sys.argv) < 2:
        print_help(commands)
    # complete
    elif sys.argv[1] == commands[0]:
        download_and_antismash(sys.argv[1])
    # download
    elif sys.argv[1] == commands[1]:
        only_download(sys.argv[1])
    # antismash
    elif sys.argv[1] == commands[2]:
        only_antismash(sys.argv[1])
    # None of the above: print help message
    else:
        print_help(commands)
        