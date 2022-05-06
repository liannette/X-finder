import os
import glob
import argparse
import subprocess
import sys
from datetime import datetime


def symbolic_links_for_antismash(outdir):
    # Create the destination directory
    if os.path.isdir(os.path.join(outdir, "antismash")) is False:
        os.mkdir(os.path.join(outdir, "antismash"))
    # Get all genomes files
    genome_files = glob.glob(os.path.join(outdir, "ncbi_dataset", "data", '*', '*.gbff'))
    # Create the symbolic links
    for genome_file in genome_files:
        new_filename = genome_file.split(os.path.sep)[-2] + '.gbff'
        dst_path = os.path.join(outdir, "antismash", new_filename)
        process = subprocess.Popen(["ln", "-s", genome_file, dst_path])


def symbolic_links_for_discoverBGC(outdir):   
    # Get all genomes files
    directories = glob.glob(os.path.join(outdir, "antismash", "*/"))
    genome_files = [x + x.split("/")[-2] + ".gbk" for x in directories]
    # Create the symbolic links 
    for genome_file in genome_files:
        dst_path = os.path.join(outdir, os.path.split(genome_file)[1])
        process = subprocess.Popen(["ln", "-s", genome_file, dst_path])


def main(mode, outdir):

    if mode == 1:
        symbolic_links_for_antismash(outdir)
    elif mode == 2:
        symbolic_links_for_discoverBGC(outdir)


if __name__ == "__main__":

    print("Commandline input:", " ".join(sys.argv))

    parser = argparse.ArgumentParser()
    parser.add_argument('--genus', action="store", dest="genus", type=str, help='Genus of the genomes')
    parser.add_argument('--mode', action="store", dest="mode", type=int, choices=[1,2], required=True, help='Specify if you want to create symbolic links for subsequent use of antismash (mode 1) of of discoverBGC (mode 2)')
    parser.add_argument('--outdir', action="store", dest="outdir", type=str, required=True, help='Directory for the symbolic links')
    
    args = parser.parse_args()
    genus = args.genus
    mode = args.mode
    outdir = args.outdir

    print("{}: Program started.".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
    main(mode=mode, outdir=outdir)
    print("{}: Program finished".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))