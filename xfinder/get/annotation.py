import os
import subprocess
import glob
import sys


def get_genomic_gbff_files(outdir):
    
    # Get all genomes files
    genomic_files = glob.glob(
        os.path.join(outdir, "ncbi_dataset", "data", '*', 'genomic.gbff')
        )
    return genomic_files
    
    # Create the symbolic links
    # renamed_files = list()
    # for genome_file in genome_files:
    #     new_filename = genome_file.split(os.path.sep)[-2] + '_genomic.gbff'
    #     dst_path = os.path.join(outdir, "antismash", new_filename)
    #     process = subprocess.Popen(["ln", "-s", genome_file, dst_path])
    #     renamed_files.append(dst_path)
    # return renamed_files


def run_antismash(outdir, genomic_files, cpus):

    # Annotate with antismash
    subprocess.call(['/data/s202633/X-finder/xfinder/get/load_antismash.sh'])
    
    # Create the destination directory
    outdir_antismash = os.path.join(outdir, "antismash")
    if os.path.isdir(outdir_antismash) is False:
        os.mkdir(outdir_antismash)
    print(f"Writing antismash results to dir: {outdir_antismash}")
    
    for gbff_file in genomic_files:
        output_basename = gbff_file.split(os.path.sep)[-2]
        print(output_basename)
        
        args = [
            "antismash", 
            "--minimal", 
            "--fullhmmer", 
            "--output-dir", outdir_antismash,
            "--output-basename", output_basename, 
            "--cpus", str(cpus), 
            gbff_file,
            ]
        output = subprocess.run(args, stdout=subprocess.PIPE, 
                                stdin=subprocess.PIPE, check=True)
        try:
            assert output.stderr == None
        except: 
            print(output.stderr)
            sys.exit(1)


def create_symbolic_links(outdir):
    # Create the symbolic links for the  annotated genomes files
    directories = glob.glob(os.path.join(outdir, "antismash", "*/"))
    genome_files = [x + x.split("/")[-2] + ".gbk" for x in directories]
    for genome_file in genome_files:
        dst_path = os.path.join(outdir, os.path.basename(genome_file))
        process = subprocess.Popen(["ln", "-s", genome_file, dst_path])


def main(outdir, cpus):
    genomic_files = get_genomic_gbff_files(outdir)
    run_antismash(outdir, genomic_files, cpus)
    create_symbolic_links(outdir)
