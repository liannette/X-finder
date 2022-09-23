import os
import sqlite3
import tempfile
import xfinder.find.mkdb
import xfinder.find.importgbk
import xfinder.find.coregen
import xfinder.find.transp
import xfinder.find.compare
import xfinder.find.cluster
from xfinder.find.common import get_git_root, get_database_dir, print_stdout, \
                                print_stderr


def make_database(database_path):
    ''' Creates the database and the tables '''

    if os.path.exists(database_path):
        os.remove(database_path)
        print_stdout("Already existing database with the same name has been " 
                     "deleted.")
    
    conn = sqlite3.connect(database_path)
    xfinder.find.mkdb.create_database_tables(conn)
    xfinder.find.mkdb.add_database_indeces(conn)
    print_stdout("New database created.")



def import_genomes_to_database(database_path, host_type, genome_dirs):
    ''' module importgbk '''

    gbk_files = xfinder.find.importgbk.list_of_gbk_files(genome_dirs)
    print_stdout("Importing {} genbank files ({}) into the database."
              "".format(len(gbk_files), host_type))

    not_imported_files = list()
    conn = sqlite3.connect(database_path)
    for i in range(len(gbk_files)):

        flag = xfinder.find.importgbk.gbk_file_to_db(conn, gbk_files[i], 
                                                     host_type)
        if flag is False:
            print_stdout("File was not imported into the database: {}"
                         "".format(gbk_files[i]))
            not_imported_files.append(gbk_files[i])

        if xfinder.find.importgbk.print_import_progress(i+1, len(gbk_files)):
            print_stdout("{}/{} files done".format(i+1, len(gbk_files)))

    conn.close()

    outfile = os.path.join(get_git_root(), "outputs",
                           "not_imported_gbk_files.txt")
    xfinder.find.importgbk.not_imported_to_file(not_imported_files, outfile)


def coregen(database_path, core_genome_path):

    print_stdout("Adding core genome information. Core genome file: {}".format(
        core_genome_path
    ))

    with tempfile.TemporaryDirectory() as temp_dir:
        
        trans_fasta = os.path.join(temp_dir, "translations.fasta")
        diamond_db = os.path.join(temp_dir, "core_genome.dmnd")
        result_file = os.path.join(temp_dir, "diamond_out.tsv")

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
    
    print_stdout("Inserted core genome information into database. {} of {} " 
                 "CDS translations aligned to the core genome with >90% "
                 "identity".format(num_core_cds, len(cds_list)))


def transp(database_path, transporter_pfams_path):

    print_stdout("Adding transporter PFAM information. Core genome "
                 "file: {}".format(transporter_pfams_path))
    transporter_pfams = xfinder.find.transp.get_transporter_pfams(
        transporter_pfams_path)
    xfinder.find.transp.add_transporter_pfam_information(
        transporter_pfams, database_path)


def compare(seed_size, gap_threshold, size_threshold, DNA_length_threshold, 
            threads, database_path):
    """ module compare """

    conn = sqlite3.connect(database_path)
    num_comparisons_start = xfinder.find.compare.count_comparisons(conn)
    host_pairs, num_query_hosts, num_ref_hosts \
        = xfinder.find.compare.get_host_comparison_pairs(conn, max_l50=3)
    conn.close()
    print_stdout("Comparing {} query hosts and {} reference hosts, which "
                 "amounts to {} comparisions.".format(
                    num_query_hosts, num_ref_hosts, len(host_pairs)))

    hit_parameters = xfinder.find.compare.HitParameters(
        seed_size, gap_threshold, size_threshold, DNA_length_threshold)

    # find hits using multithreading
    xfinder.find.compare.find_hits_multithread(host_pairs, hit_parameters, 
                                               database_path, threads)

    # Calculate how many new comparions were made 
    conn = sqlite3.connect(database_path)
    comparison_count_new = xfinder.find.compare.count_comparisons(conn) \
                           - num_comparisons_start
    conn.close()
    print_stdout("{} new comparisons have been added to "
        "the database.".format(comparison_count_new))


def cluster(threads, database_path):
    print_stdout("Grouping similar hits into clusters")
    conn = sqlite3.connect(database_path)
    clustered_hits = xfinder.find.cluster.cluster_hits(conn, threads, 
                                                       database_path)
    xfinder.find.cluster.add_cluster_information_to_db(clustered_hits, conn)
    conn.close()

#############################################################################
def main():
    database = "test_db.db"
    database_path = os.path.join(get_database_dir(), database)

    # Create db
    make_database(database_path)

    # Import genomes
    host_type = "ref"
    genome_dirs = [os.path.join(get_git_root(), "data", "genomes", "ref_genomes")]
    import_genomes_to_database(database_path, host_type, genome_dirs)
    host_type = "query"
    genome_dirs = [os.path.join(get_git_root(), "data", "genomes", "query_genomes")]
    import_genomes_to_database(database_path, host_type, genome_dirs)

    # Add core genome information
    core_genome_path = os.path.join(get_git_root(), "data", "core genome", 
                                    "Streptomyces.fasta") 
    coregen(database_path, core_genome_path)

    # Add transporter pfam information
    file_path = os.path.join(get_git_root(), "data", "transporter pfams",
                             "all_transporter_pfams.txt")
    transp(database_path, file_path)

    # Run comparison to find hits
    seed_size = 2
    gap_threshold = 2
    size_threshold = 6
    DNA_length_threshold = 7000
    threads = 2
    compare(seed_size, gap_threshold, size_threshold, DNA_length_threshold, 
                threads, database_path)

    # cluster hits
    cluster(threads, database_path)


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

"""


