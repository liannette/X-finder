"""
xfinder

Find novel biosynthetic gene cluster
"""

__version__ = "0.1.0"
__author__ = 'Annette Lien'

from xfinder.mkdb import make_database
from xfinder.importgbk import import_genomes
from xfinder.coregen import add_coregenome_info
from xfinder.transp import add_transporter_info
from xfinder.compare import compare_genomes
from xfinder.cluster import cluster_hits
from xfinder.results import export_results


def run_all(database_path, ref_genome_dirs, query_genome_dirs, 
            core_genome_path, transporter_pfams_path, seed_size, 
            gap_threshold, size_threshold, DNA_length_threshold, max_l50,
            threads, core_genome_cutoff, transporter_cutoff, out_dir):
        
    # Create db
    make_database(database_path, out_dir)

    # Import genomes
    # I need to add a check that the directories indeed exist
    import_genomes(database_path, "ref", ref_genome_dirs, out_dir)
    import_genomes(database_path, "query", query_genome_dirs, out_dir)

    # Add core genome and transporter pfam information
    add_coregenome_info(database_path, core_genome_path, out_dir)
    add_transporter_info(database_path, transporter_pfams_path, out_dir)

    # find hits
    compare_genomes(seed_size, gap_threshold, size_threshold, 
                    DNA_length_threshold, threads, max_l50, database_path, 
                    out_dir)
    # cluster hits
    cluster_hits(threads, database_path, out_dir)

    # print results
    export_results(core_genome_cutoff, transporter_cutoff, database_path, 
                    out_dir)
    