import os
import xfinder 

out_dir = "/data/s202633/X-finder_outputs/myxococcales_complete"


database_path = os.path.join(out_dir, "database.db")

genomes_dir = "/data/s202633/X-finder_inputs/genomes"
ref_genome_dirs = [
    "/data/s202633/X-finder_inputs/genomes/streptomyces_internal", 
    "/data/s202633/X-finder_inputs/genomes/actinobacteria_internal",
    ]    
query_genome_dirs = [
    "/data/s202633/X-finder_inputs/genomes/myxococcales_complete",
    ]    

# out_dir = "/data/s202633/X-finder_outputs/test"
# database = "test_db.db"
# database_path = os.path.join(out_dir, database)
# genomes_dir = "/data/s202633/X-finder_inputs/genomes"
# ref_genome_dirs = [os.path.join(genomes_dir, "ref_genomes"),]
# query_genome_dirs = [os.path.join(genomes_dir, "query_genomes"),]  


core_genome_path = "/data/s202633/X-finder_inputs/core_genome/Streptomyces.fasta"
transporter_pfams_path = "/data/s202633/X-finder_inputs/transporter_pfams/all_transporter_pfams.txt"

seed_size = 2
gap_threshold = 2
size_threshold = 6
DNA_length_threshold = 7000
threads = 8
max_l50 = 3

core_genome_cutoff = 0.5
transporter_cutoff = 0.2


xfinder.run_all(
    database_path, ref_genome_dirs, query_genome_dirs, 
    core_genome_path, transporter_pfams_path, seed_size, 
    gap_threshold, size_threshold, DNA_length_threshold, max_l50,
    threads, core_genome_cutoff, transporter_cutoff, out_dir
    )
