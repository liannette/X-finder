import os
from subprocess import Popen, PIPE
import glob
from getgenomes.common import print_log, print_stdout_stderr, create_dir
from pathlib import Path

def _run_antismash(genome_file, outdir, outdir_antismash, threads, 
                   antismash_path):
    
    args = [str(antismash_path), 
            "--minimal", "--fullhmmer", 
            "--cpus", str(threads), 
            "--output-dir", 
            str(os.path.join(outdir_antismash, Path(genome_file).name)), 
            str(genome_file)]
    process = Popen(args, stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate() 
    if process.returncode != 0: 
        print_log("An error occured while running the following command: "
                  f"{' '.join(args)}", outdir)
    print_stdout_stderr(stdout, stderr, outdir)


def _create_links(outdir_antismash, outdir):
    """ Create the symbolic links for the annotated genomes files """
    
    # Destination directory for the xfinder input genomes
    dst_dir = os.path.join(outdir_antismash, "all")
    if os.path.isdir(dst_dir) is False:
        os.mkdir(dst_dir)
    
    for gbk_file in glob.glob(os.path.join(outdir_antismash, "*", "*.gbk")):
        # no region gbk files 
        if Path(gbk_file).stem in Path(gbk_file).parent.name:
            # Create the link
            process = Popen(["ln", gbk_file, dst_dir], stdout=PIPE, 
                            stderr=PIPE)
            stdout, stderr = process.communicate() 
            print_stdout_stderr(stdout, stderr, outdir)


def run_antismash(indir, outdir, threads, antismash_path):
    """ Runs antiSMASH on all file downloaded previouly with 
    ncbi-genome-download. Settings for antiSMASH are --minimal and 
    --fullhmmer. Creates a directory with the name 'final' in the output
    directory, containing symbolic links to the genome genbank files 
    from the antiSMASH output."""

    # Check the arguments
    if os.path.isdir(indir) is False:
        raise RuntimeError(f"Not a directory: {indir}")
    create_dir(outdir)

    # Create the destination directory
    outdir_antismash = os.path.join(outdir, "antismash")
    create_dir(outdir_antismash)

    # Get a list of all paths of the genome files
    genome_files = glob.glob(os.path.join(indir, "*"))
    print_log(f"Started antiSMASH on {len(genome_files)} files", outdir)
    for genome_file in genome_files:
        _run_antismash(genome_file, outdir, outdir_antismash, threads, 
                       antismash_path)
        
    # Create a directory with symbolic links of all antismash gbk genomes   
    _create_links(outdir_antismash, outdir)
    