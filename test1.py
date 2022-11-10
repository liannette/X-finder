import os
import xfinder 

# out_dir = "/data/s202633/X-finder_outputs/myxococcales_complete"
# ref_genome_dirs = [
#     "/data/s202633/X-finder_inputs/genomes/streptomyces_internal", 
#     "/data/s202633/X-finder_inputs/genomes/actinobacteria_internal",
#     ]    
# query_genome_dirs = [
#     "/data/s202633/X-finder_inputs/genomes/myxococcales_complete",
#     ]    

out_dir = "/data/s202633/X-finder_outputs/test"
genomes_dir = "/data/s202633/X-finder_inputs/genomes"
ref_genome_dirs = ["/data/s202633/X-finder_inputs/genomes/ref_genomes"]
query_genome_dirs = ["/data/s202633/X-finder_inputs/genomes/query_genomes"]



database_path = os.path.join(out_dir, "database.db")
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


# xfinder.run_all(
#     database_path, ref_genome_dirs, query_genome_dirs, 
#     core_genome_path, transporter_pfams_path, seed_size, 
#     gap_threshold, size_threshold, DNA_length_threshold, max_l50,
#     threads, core_genome_cutoff, transporter_cutoff, out_dir
#     )



# # Create db
# xfinder.make_database(database_path, out_dir)

# # Import genomes
# # I need to add a check that the directories indeed exist
# xfinder.import_genomes(database_path, "ref", ref_genome_dirs, out_dir)
# xfinder.import_genomes(database_path, "query", query_genome_dirs, out_dir)

# Add core genome and transporter pfam information
xfinder.add_coregenome_info(database_path, core_genome_path, out_dir)
xfinder.add_transporter_info(database_path, transporter_pfams_path, out_dir)

# find hits
xfinder.compare_genomes(seed_size, gap_threshold, size_threshold, 
                DNA_length_threshold, threads, max_l50, database_path, 
                out_dir)
# cluster hits
xfinder.cluster_hits(threads, database_path, out_dir)

# print results
xfinder.export_results(core_genome_cutoff, transporter_cutoff, database_path, 
                out_dir)