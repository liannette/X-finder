import os
import sqlite3
import tempfile
import traceback
import xfinder.find.mkdb
import xfinder.find.importgbk
import xfinder.find.coregen
import xfinder.find.transp
import xfinder.find.compare
import xfinder.find.cluster
import xfinder.find.results
from xfinder.find.common import get_git_root, print_stdout, print_stderr


def make_database(database_path, out_dir):
    ''' Creates the database and the tables '''

    if os.path.exists(database_path):
        os.remove(database_path)
        print_stdout("Already existing database with the same name has been " 
                     "deleted.", out_dir)
    
    conn = sqlite3.connect(database_path)
    xfinder.find.mkdb.create_database_tables(conn)
    xfinder.find.mkdb.add_database_indeces(conn)
    print_stdout("New database created.", out_dir)



def import_genomes_to_database(database_path, host_type, genome_dirs, out_dir):
    ''' module importgbk '''

    gbk_files = xfinder.find.importgbk.list_of_gbk_files(genome_dirs)
    print_stdout(f"Importing {len(gbk_files)} genbank files ({host_type}) "
                 "into the database.", out_dir)
    
    not_imported_files = list()
    conn = sqlite3.connect(database_path)
    for i in range(len(gbk_files)):
        success = xfinder.find.importgbk.gbk_file_to_db(conn, gbk_files[i], 
                                                     host_type)
        if success is False:
            print_stdout("File was not imported into the database: "
                         f"{gbk_files[i]}", out_dir)
            not_imported_files.append(gbk_files[i])
        if xfinder.find.importgbk.print_import_progress(i+1, len(gbk_files)):
            print_stdout("{}/{} files done".format(i+1, len(gbk_files)), 
                         out_dir)
    conn.close()

    xfinder.find.importgbk.not_imported_to_file(
        not_imported_files, 
        outfile_path=os.path.join(out_dir, "not_imported_gbk_files.txt"))


def coregen(database_path, core_genome_path, out_dir):

    print_stdout("Adding core genome information. Core genome file: "
                 f"{core_genome_path}", out_dir)

    with tempfile.TemporaryDirectory() as temp_dir:
        
        trans_fasta = os.path.join(temp_dir, "translations.fasta")
        diamond_db = os.path.join(temp_dir, "core_genome.dmnd")
        result_file = os.path.join(temp_dir, "diamond_out.tsv")

        #Rewrite this so das cds_list is not nesseccary!

        # Write a fasta file with all cds translations
        cds_list = xfinder.find.coregen.cds_translations_to_fasta(
            database_path, trans_fasta)
        # Create diamond databank of the core genome
        xfinder.find.coregen.create_diamond_database(
            diamond_db, core_genome_path)
        # Run diamond
        xfinder.find.coregen.run_diamond(trans_fasta, diamond_db, result_file)

        # Add the core genome information to the database
        num_core_cds = xfinder.find.coregen.add_core_genome_information(
            database_path, result_file, cds_list)
    
    print_stdout("Inserted core genome information into database. "
                 f"{num_core_cds}/{len(cds_list)} CDS translations aligned " 
                 "to the core genome with >90% identity", out_dir)


def transp(database_path, transporter_pfams_path, out_dir):

    print_stdout("Adding transporter PFAM information. Core genome file: "
                 f"{transporter_pfams_path}", out_dir)
    transporter_pfams = xfinder.find.transp.get_transporter_pfams(
        transporter_pfams_path)
    xfinder.find.transp.add_transporter_pfam_information(
        transporter_pfams, database_path)


def compare(seed_size, gap_threshold, size_threshold, DNA_length_threshold, 
            threads, max_l50, database_path, out_dir):
    """ module compare """

    conn = sqlite3.connect(database_path)
    
    # Get the number of comparisons already in the database
    num_comparisons_start = xfinder.find.compare.count_comparisons(conn)
    
    # Get the hostIDs
    query_hostIDs = xfinder.find.compare.get_host_ids(
        conn, 
        "query", 
        max_l50
        )
    ref_hostIDs = xfinder.find.compare.get_host_ids(
        conn, 
        "ref", 
        max_l50
        )
        
    conn.close()
    
    print_stdout(
        f"Comparing {len(query_hostIDs)} query hosts and {len(ref_hostIDs)} "
        "reference hosts, which amounts to  "
        f"{len(query_hostIDs)*len(ref_hostIDs)} comparisions.", out_dir)


    # Put the hit parameter in a nice class
    hit_parameters = xfinder.find.compare.HitParameters(
        seed_size, 
        gap_threshold, 
        size_threshold, 
        DNA_length_threshold
        )

    # find hits using multithreading
    xfinder.find.compare.find_hits_multithread(
        query_hostIDs, 
        ref_hostIDs,
        hit_parameters, 
        database_path, 
        threads
        )

    # Calculate how many new comparions were made 
    conn = sqlite3.connect(database_path)
    comparison_count_new = xfinder.find.compare.count_comparisons(conn) \
                           - num_comparisons_start
    conn.close()
    print_stdout(f"{comparison_count_new} new comparisons have been added to "
                 "the database.", out_dir)


def cluster(threads, database_path, out_dir):
    print_stdout("Grouping similar hits into clusters", out_dir)
    clustered_sublists = xfinder.find.cluster.cluster_hits(
        threads, 
        database_path
        )
    print_stdout(f"{len(clustered_sublists)} cluster found in total", out_dir)
    xfinder.find.cluster.add_cluster_to_db(
        clustered_sublists, 
        database_path
        )


def results(core_genome_cutoff, transporter_cutoff, database_path, out_dir):
    
    cluster_list = xfinder.find.results.get_filtered_cluster(
        core_genome_cutoff, 
        transporter_cutoff, 
        database_path
        )
    
    print_stdout(f"Printing results. {len(cluster_list)} cluster fullfill the "
                 f"cutoff criteria (core genome: {core_genome_cutoff}, "
                 f"transporter: {transporter_cutoff})", out_dir)
                 
    xfinder.find.results.write_result_files(
        cluster_list, 
        database_path, 
        core_genome_cutoff,
        transporter_cutoff, 
        out_dir
        )


def run_xfinder(database_path, ref_genome_dirs, query_genome_dirs, 
                core_genome_path, transporter_pfams_path, seed_size, 
                gap_threshold, size_threshold, DNA_length_threshold, max_l50,
                threads, core_genome_cutoff, transporter_cutoff, out_dir):
    

        
    # # Create db
    # make_database(database_path, out_dir)

    # # Import genomes
    # # I need to add a check that the directories indeed exist
    # import_genomes_to_database(database_path, "ref", ref_genome_dirs, out_dir)
    import_genomes_to_database(database_path, "query", query_genome_dirs, out_dir)

    # Add core genome and transporter pfam information
    coregen(database_path, core_genome_path, out_dir)
    transp(database_path, transporter_pfams_path, out_dir)

    # find hits
    compare(seed_size, gap_threshold, size_threshold, DNA_length_threshold, 
                threads, max_l50, database_path, out_dir)
    # cluster hits
    cluster(threads, database_path, out_dir)

    # print results
    results(core_genome_cutoff, transporter_cutoff, database_path, out_dir)



#############################################################################


def main():
    
    out_dir = "/data/s202633/X-finder_outputs/myxococcales_complete"
    database = "myxococcales_complete.db"
    
    database_path = os.path.join(out_dir, database)
    
    genomes_dir = "/data/s202633/X-finder_inputs/genomes"
    ref_genome_dirs = [
        os.path.join(genomes_dir, "streptomyces_internal"), 
        os.path.join(genomes_dir, "actinobacteria_internal"), 
        ]    
    query_genome_dirs = [
        os.path.join(genomes_dir, "myxococcales_complete"),
        ]    

    # out_dir = "/data/s202633/X-finder_outputs/test"
    # database = "test_db.db"
    # database_path = os.path.join(out_dir, database)
    # genomes_dir = "/data/s202633/X-finder_inputs/genomes"
    # ref_genome_dirs = [os.path.join(genomes_dir, "ref_genomes"),]
    # query_genome_dirs = [os.path.join(genomes_dir, "query_genomes"),]  


    core_genome_path = os.path.join(get_git_root(), "data", "core genome", 
                                    "Streptomyces.fasta") 

    transporter_pfams_path = os.path.join(get_git_root(), "data",
                             "transporter pfams", "all_transporter_pfams.txt")

    seed_size = 2
    gap_threshold = 2
    size_threshold = 6
    DNA_length_threshold = 7000
    threads = 8
    max_l50 = 3
    
    core_genome_cutoff = 0.5
    transporter_cutoff = 0.2


    try:
        # Create output directory
        # if os.path.isdir(out_dir):
        #     out_dir = None
        #     raise RuntimeError("Output directory already exists, aborting for "
        #                        "safety")
        # else:
        #     os.mkdir(out_dir)
            
        # run xfinder
        run_xfinder(
            database_path, ref_genome_dirs, query_genome_dirs, 
            core_genome_path, transporter_pfams_path, seed_size, 
            gap_threshold, size_threshold, DNA_length_threshold, max_l50,
            threads, core_genome_cutoff, transporter_cutoff, out_dir)
    except:
        print_stderr(traceback.format_exc(), out_dir)


if __name__ == "__main__":
    main()


"""
## importgbk
if __name__ == "__main__":

    print("Commandline input:", " ".join(sys.argv))

    parser = argparse.ArgumentParser()
    parser.add_argument('-db', action="store", dest="database", type=str, default="database.db", help='Name of the sqlite database (default=database.db')
    parser.add_argument('--test', action="store_true", help="creates a test database with only 4 hosts")
    parser.add_argument('-q', nargs="+", action="store", dest="query_dirs", default="", help='Directory name of the query host genomes, multiple directories possible seperated with spaces. Directory must be located in discoverGBC/data/genomes and genomes must be in genbank format')
    parser.add_argument('-r', nargs="+", action="store", dest="ref_dirs", default="", help='Directory name of the ref host genomes, multiple directories possible seperated with spaces. Directory must be located in discoverGBC/data/genomes and genomes must be in genbank format')

    args = parser.parse_args()
    database = args.database
    test_flag = args.test
    query_dirs = args.query_dirs
    ref_dirs = args.ref_dirs

    print("{}: Program started. Database: {}".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), database))
    main(test_flag=test_flag, database=database, query_dirs=query_dirs, ref_dirs=ref_dirs)
    print("{}: Program finished".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))

## coregen
if __name__ == "__main__":

    print("Commandline input:", " ".join(sys.argv))

    parser = argparse.ArgumentParser()
    parser.add_argument('-db', action="store", dest="database", type=str, default="database.db", help='Name of the sql database (default=database.db)')
    args = parser.parse_args()
    database = args.database

    print("{}: Program started. Database: {}".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), database))
    main(database=database)
    print("{}: Program finished".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))

## compare
if __name__ == "__main__":

    print("Commandline input:", " ".join(sys.argv))

    parser = argparse.ArgumentParser()
    parser.add_argument('-t', action="store", dest="threads", type=int, default=1, help='Number of threads for multiprocessing (default=1)')
    parser.add_argument('-db', action="store", dest="database", type=str, default="database.db", help='Name of the sql database (default=database.db')
    args = parser.parse_args()
    threads = args.threads
    database = args.database

    print("{}: Program started. Threads: {}, database: {}".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), threads, database))
    main(seed_size=2, 
         gap_threshold=2, 
         size_threshold=6, 
         DNA_length_threshold=7000, 
         threads=threads, 
         database=database
         )
    print("{}: Program finished".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
    # to check if the mode is the main program. if it is the if condition is true then the main() will be excuted
    # if it is mode imported by other mode or the main program. __name__= the name of this file.
    #the if condition will be faulse, and the main() will not be excuted
    # https://stackoverflow.com/questions/419163/what-does-if-name-main-do

    # Regarding multithreading: SQLite accepts several connections at the same time, but only one connection can write at the same time.

## cluster1
if __name__ == "__main__":

    print("Commandline input:", " ".join(sys.argv))

    parser = argparse.ArgumentParser()
    parser.add_argument('-t', action="store", dest="threads", type=int, default=1, help='Number of threads for multiprocessing (default=1)')
    parser.add_argument('-db', action="store", dest="database", type=str, default="database.db", help='Name of the sql database (default=database.db)')

    args = parser.parse_args()
    threads = args.threads
    database = args.database

    print("{}: Program started. Database: {}, Threads: {}".format(
        datetime.now().strftime("%d/%m/%Y %H:%M:%S"), database, threads))
    main(threads=threads, database=database)
    print("{}: Program finished".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))


# results =================================================================

if __name__ == "__main__":

    print("Commandline input:", " ".join(sys.argv))

    parser = argparse.ArgumentParser()
    parser.add_argument('-db', action="store", dest="database", type=str, default="database.db", help='Name of the sql database (default=database.db')
    parser.add_argument('-core', action="store", dest="core_genome_cutoff", type=float, default=0.5, help='Cutoff fot the core genome indicator. From all ref sublists in the hits, is the max fraction of cds that aligned to the core genome. (default=0.5)')
    parser.add_argument('-transp', action="store", dest="transporter_cutoff", type=float, default=0.2, help='Cutoff for the transporter indicatior. Value indicates the fraction of transporter pfams in the core pfams of the cluster. (default=0.2)')
    parser.add_argument('-o', action="store", dest="outfile_prefix", type=str, default="results", help='Outfile suffix for the summary and detailed results. (default=results)')
    
    args = parser.parse_args()
    database = args.database
    core_genome_cutoff = args.core_genome_cutoff
    transporter_cutoff = args.transporter_cutoff
    outfile_prefix = args.outfile_prefix

    print("{}: Program started. Database: {}, core_genome_cutoff: {}, transporter_cutoff: {}, outfile_prefix: {}" \
        .format(
            datetime.now().strftime("%d/%m/%Y %H:%M:%S"), 
            database,
            core_genome_cutoff,
            transporter_cutoff,
            outfile_prefix
            )
        )
    main(database=database, core_genome_cutoff=core_genome_cutoff, transporter_cutoff=transporter_cutoff, outfile_prefix=outfile_prefix)
    print("{}: Program finished".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))

"""


