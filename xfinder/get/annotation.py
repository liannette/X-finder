import os
import subprocess
import glob
import sys

def symbolic_links_for_antismash(outdir):
    # Create the destination directory
    if os.path.isdir(os.path.join(outdir, "antismash")) is False:
        os.mkdir(os.path.join(outdir, "antismash"))
    # Get all genomes files
    genome_files = glob.glob(os.path.join(outdir, "ncbi_dataset", "data", '*', 'genomic.gbff'))
    # Create the symbolic links
    renamed_files = list()
    for genome_file in genome_files:
        new_filename = genome_file.split(os.path.sep)[-2] + '_genomic.gbff'
        dst_path = os.path.join(outdir, "antismash", new_filename)
        process = subprocess.Popen(["ln", "-s", genome_file, dst_path])
        renamed_files.append(dst_path)
    return renamed_files


def run_antismash(outdir, genomic_files, cpus):

    # Annotate with antismash
    subprocess.call(['./load_antismash.sh'])
    for gbk_file in genomic_files:
        args = [
            "antismash", "--minimal", "--fullhmmer", "--cpus", str(cpus), gbk_file
            ]
        output = subprocess.run(args, stdout=subprocess.PIPE, stdin=subprocess.PIPE, check=True)
        try:
            assert output.stderr == None
        except: 
            print(output.stderr)
            sys.exit(1)

    # Create the symbolic links for the  annotated genomes files
    directories = glob.glob(os.path.join(outdir, "antismash", "*/"))
    genome_files = [x + x.split("/")[-2] + ".gbk" for x in directories]
    for genome_file in genome_files:
        dst_path = os.path.join(outdir, os.path.split(genome_file)[1])
        process = subprocess.Popen(["ln", "-s", genome_file, dst_path])


def main(outdir, cpus):
    genomic_files = symbolic_links_for_antismash(outdir)
    run_antismash(outdir, genomic_files, cpus)
