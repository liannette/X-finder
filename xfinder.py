"""
Author: Annette Lien (a.lien@posteo.de), Xinglin Jiang

usage:
python xfinder.py -h
"""

import argparse
import xfinder
import os
from multiprocessing import cpu_count
import glob


def get_commands():
    parser = argparse.ArgumentParser(
        description="")
    parser.add_argument(
        "-o", "--output_dir", dest="out_dir", required=True,
        help="Output directory, this will contain all output data files. "
        "Must not exist already.", metavar="<dir>")
    parser.add_argument(
        "-q", "--queries_dirs", dest="query_genome_dirs", help="Input "
        "directories of query gbk files, Genomes must be in Genbank format "
        "with '.gbk' file extension and contain antismash annotation "
        "(--fullhmmer). Multiple directories are possible, paths must be "
        "seperated with a ','", required=True, metavar="<dirs>")
    parser.add_argument(
        "-r", "--ref_dirs", dest="ref_genome_dirs", help="Input "
        "directories of reference genomes. Genomes must be in Genbank format "
        "with '.gbk' file extension and contain antismash annotation "
        "(--fullhmmer). Multiple directories are possible, paths must be "
        "seperated with a ','", required=True, metavar="<dirs>")
    parser.add_argument(
        "--core_genome", dest="core_genome_path", required=True, 
        metavar="<file>", help="Fasta file containing core genome sequences.")
    parser.add_argument(
        "--transporter_pfams", dest="transporter_pfams_path", default=None,
        metavar="<file>", help=" File containing all PFAM domains that are "
        " associated with transporter function, one PFAM number (without the "
        "PF) per line. (default: None)")
    parser.add_argument(
        "-t", "--threads", dest="threads", default=cpu_count(),
        help="Set the number of threads the script may use (default: use all \
        available cores)", type=int, metavar="<int>")
    parser.add_argument(
        "--max_l50", dest="max_l50", default=None, help="Only considers input "
        "genomes with a L50 value of this or lower. L50 is defined as count "
        "of smallest number of contigs whose length sum makes up half of "
        "genome size (default: None)", type=int, metavar="<int>")
    parser.add_argument(
        "--min_seed_size", dest="seed_size", default=2, help="The minimum "
        "number of pfam domains that match exactly during the initial step of "
        "finding common pfam domains between a query and ref genome "
        "(each genome represented by a genomepfamsequence). Do not change "
        "this value, except you know what you are doing! (default: 2)",
        type=int, metavar="<int>")
    parser.add_argument(
        "--max_gap", dest="gap_threshold", default=2, help="The maximum number"
        " of pfam domains in a gap, when combining neighboring/overlapping "
        "exact matches (default: 2)", type=int, metavar="<int>")
    parser.add_argument(
        "--min_sublist_size", dest="size_threshold", default=6, help="Both "
        "sublists of a hit must contain at least this many PFAM domains. "
        "(default: 6)", type=int, metavar="<int>")
    parser.add_argument(
        "--min_dna_length", dest="DNA_length_threshold", default=7000, 
        help="The DNA sequence length of the query sublist must be at this "
        "long (default: 7000)", type=int, metavar="<int>")
    parser.add_argument(
        "--core_cutoff", dest="core_genome_cutoff", default=0.5, help="Cutoff "
        "for core genome filtering (default: 0.5)", type=float, 
        metavar="<float>")
    parser.add_argument(
        "--transporter_cutoff", dest="transporter_cutoff", default=0.2, 
        help="Cutoff for transporter function filtering (default: 0.2)", 
        type=float, metavar="<float>")
    return parser.parse_args()


def check_if_files_exist(out_dir, ref_genome_dirs, query_genome_dirs, 
                         core_genome_path, transporter_pfams_path):
    if os.path.isdir(out_dir):
        raise RuntimeError("Output directory already exists, aborting for "
                           "safety")
    for genome_dir in ref_genome_dirs:
        if os.path.exists(genome_dir) is False:
            raise RuntimeError("Aborting because input dir for reference gbk "
                               f"files does not exist: {genome_dir}")
        if len(glob.glob(os.path.join(genome_dir, "*.gbk"))) == 0:
            raise RuntimeError("Aborting because input dir for reference gbk "
                               f"files contains no gbk files: {genome_dir}")
    for genome_dir in query_genome_dirs:
        if os.path.exists(genome_dir) is False:
            raise RuntimeError("Aborting because input dir for query gbk "
                               f"files does not exist: {genome_dir}")
        if len(glob.glob(os.path.join(genome_dir, "*.gbk"))) == 0:
            raise RuntimeError("Aborting because input dir for query gbk "
                               f"files contains no gbk files: {genome_dir}")
    if os.path.exists(core_genome_path) is False:
        raise RuntimeError("Aborting because core genome file does not exist: "
                           f"{core_genome_path}")
    if transporter_pfams_path is not None:    
        if os.path.exists(transporter_pfams_path) is False:
            raise RuntimeError("Aborting because transporter pfams file does "
                           f"not exist: {transporter_pfams_path}")
    

if __name__ == "__main__":
    
    args = get_commands()
    
    check_if_files_exist(
        args.out_dir, 
        args.ref_genome_dirs.split(","), 
        args.query_genome_dirs.split(","),
        args.core_genome_path,
        args.transporter_pfams_path,
        )
    
    xfinder.run_all(
        database_path=os.path.join(args.out_dir, "database.db"), 
        ref_genome_dirs=args.ref_genome_dirs.split(","), 
        query_genome_dirs=args.query_genome_dirs.split(","), 
        core_genome_path=args.core_genome_path, 
        transporter_pfams_path=args.transporter_pfams_path, 
        seed_size=args.seed_size, 
        gap_threshold=args.gap_threshold, 
        size_threshold=args.size_threshold, 
        DNA_length_threshold=args.DNA_length_threshold, 
        max_l50=args.max_l50, 
        threads=args.threads, 
        core_genome_cutoff=args.core_genome_cutoff, 
        transporter_cutoff=args.transporter_cutoff, 
        out_dir=args.out_dir, 
        )
    
    
    