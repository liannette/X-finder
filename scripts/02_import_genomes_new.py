# Reference http://www2.warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/genbank/
# Reference http://blog.csdn.net/y2701310012/article/details/42218255
# Here we extract all PFAM numbers and corresponding locations and locus_tags from an annotated genbank file.
# For pfam features that overlap with eachother, only the one with smallest pfumnumber will be kept. (alternatively we can choose the one with highest score (not done yet))
# however short overlap smaller than 1/2 of both pfum features will be tolerated.
# we did not consider the drection of pfam, because we want to tolerate gene inversion.

import sys
import argparse
import glob
import sqlite3
from Bio import SeqIO
from general_functions import connect_to_db, get_base_dir
import os.path
import re 
import traceback
from datetime import datetime


# Functions ====================================================================


def get_hosttype(host, ref_genus):
    """
    Checks if the host starts with the ref_genus
    """
    if host.startswith(ref_genus):
        hosttype = "reference"
    else:
        hosttype = "query" 
    return hosttype


def insert_host_to_DB(conn, host, description):
    '''
    insert a new host into the hosts table
    :param conn:
    :param host:
    :param description:
    '''
    c = conn.cursor()
    sql =   ''' 
            INSERT OR IGNORE INTO hosts (host, description, hosttype) 
                           VALUES (?,?,?) 
            '''
    c.execute(sql, (host, description, get_hosttype(description, "Streptomyces")))
    c.close()


def insert_contig_to_DB(conn, host, accession, sequence_length):
    '''
    insert a new contig into the contigs table
    '''
    c = conn.cursor()
    sql =   ''' 
            INSERT OR IGNORE INTO contigs (contig, host, sequence_length) 
                           VALUES (?,?,?) 
            '''
    c.execute(sql, (accession, host, sequence_length))
    c.close()


def insert_cds_to_DB(conn, cds):
    '''
    insert the cds locus_tag, product and translation into the table cds
    :param conn:
    :param cds_products:
    :return:    
    '''
    c = conn.cursor()
    sql = ''' INSERT INTO cds (contig, locus_tag, product, translation, cds_start, cds_end)
              VALUES (?,?,?,?,?,?) '''

    c.executemany(sql, cds)
    c.close()


def insert_gpfamsequence_to_DB(conn, gpfamsequence):
    '''
    insert a new gpfam into the pfam table
    gpfams from a genomes make up a gpfamsequence
    :param conn:
    :param gpfamsequence:
    :return:
    '''
    c = conn.cursor()
    sql = ''' INSERT INTO pfams (locus_tag, pfamnumber, pfamstart, pfamend, strand, contig)
              VALUES(?,?,?,?,?,?) '''
    c.executemany(sql, gpfamsequence)
    c.close()


def extract_from_gbfile_to_DB(genbank_file, host, conn):
    '''
    host(genome sequenceID, etc) and contigs are extracted and inserted to DB; 
    pfams [(pfamnumber, start, end, locus_tag)] are extracted from the sequence and make up a gpfamsequence and is then inserted to DB; 
    if two pfams overlap with eachother, and the overlap area is bigger then half of either of them, only the one with smallest pfumnumber will be kept.
    '''
    for rec in SeqIO.parse(genbank_file, "genbank"):
        # Some PFAM_domain features don't have PFAM ID (this is actually a  bug from antismash (they use a dictionary of pfam name and pfam ID which need to be updated mannually), remember to check if this happen too often, if yes ask Kai to fix it).so we used try expect. .
        description = rec.description
        contig_accession = rec.id
        insert_host_to_DB(conn, host, description)
        insert_contig_to_DB(conn, host, contig_accession, len(rec.seq))

        cds_list = []
        gpfamsequence = []
        pfamstart = 0
        pfamend = 0
        pfamnumber = 10000     #just a number bigger than any real pfamnumber.
        if rec.features:
            for feature in rec.features:
                if feature.type == "CDS":
                    try:
                        locus_tag = feature.qualifiers["locus_tag"][0]
                        cds_product = feature.qualifiers["product"][0]
                        translation = feature.qualifiers["translation"][0]
                        cdsstart = feature.location.nofuzzy_start
                        cdsend = feature.location.nofuzzy_end
                        cds_list.append( (contig_accession, locus_tag, cds_product, translation, cdsstart, cdsend) )
                    except:
                        continue
                elif feature.type == "PFAM_domain":
                    try:
                        apfamid = feature.qualifiers["db_xref"][0]
                        apfamnumber = int(re.search(r"PF(\d{5})\D", apfamid).group(1))
                        apfam_locus_tag = feature.qualifiers["locus_tag"][0]
                        apfamsize = feature.location.nofuzzy_end - feature.location.nofuzzy_start
                        astrand = feature.strand
                        overlap = pfamend - feature.location.nofuzzy_start  #nofuzzy_start left most (minimum) value, regardless of strand
                    except:
                        continue #pfamnumber can not be found
                    if overlap > apfamsize/2 or overlap > (pfamend - pfamstart)/2:
                        if apfamnumber >= pfamnumber:
                            #print ("overlap found")
                            continue

                        #remove the previous pfam
                        gpfamsequence.pop()
                        #print ("overlapped pfam found and removed")

                    pfamnumber = apfamnumber
                    pfam_locus_tag = apfam_locus_tag
                    pfamstart = feature.location.nofuzzy_start
                    pfamend = feature.location.nofuzzy_end
                    strand = astrand
                    pfam = (pfam_locus_tag, pfamnumber, pfamstart, pfamend, strand, contig_accession)
                    gpfamsequence.append(pfam)

        #print ("in total ", len(cds_list), " CDS found")
        #print ("in total ", len(gpfamsequence), " pfams found")

        insert_cds_to_DB(conn, cds_list)
        insert_gpfamsequence_to_DB(conn, gpfamsequence)


def import_gbk_files_to_DB(database, test_flag):
    '''
    '''
    # Get all genome files
    genomes_dir = os.path.join(get_base_dir(), "data", "genomes", "Actinobacteria_internal", "*.gbk")
    genbank_files = glob.glob(genomes_dir)
    for genbank_file in genbank_files:
        print(genbank_file)
    sys.exit(1)
    if test_flag:
        genbank_files = genbank_files[:20]
    print("{}: Importing {} genomes into the database '{}'.".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), len(genbank_files), database))

    # Extract information from the genbank files into the database
    conn = connect_to_db(database)
    missing_files = list()
    for i in range(len(genbank_files)):
        genbank_file = genbank_files[i]
        try:
            extract_from_gbfile_to_DB(genbank_file, hostID, conn)
        except:
            missing_files.append(genbank_file)
            print("{}: Following file could not be imported to the database: {}".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), genbank_file))
            traceback.print_exc()
        if i+1 % 20 == 0 or i+1 == len(directories):
            print("{}: {}/{} genomes done".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), i+1, len(directories)))
    conn.close()

    # print list of genomes that were not imported
    if test_flag:
        outfile = open(os.path.join(get_base_dir(), "results", "not_imported_genome_files_test.txt"), "w")
    else: 
        outfile = open(os.path.join(get_base_dir(), "results", "not_imported_genome_files.txt"), "w")
    print("Following files were not imported into the database:", file=outfile)
    for filename in missing_files:
        print(filename, file=outfile)
    outfile.close()


def read_hosts_from_db(conn):
    '''
    reads all hosts from table hosts
    :param conn:
    :return: a python list of hosts
    '''
    c = conn.cursor()
    rows = c.execute(''' SELECT host FROM hosts''',).fetchall()
    c.close()
    return rows


def read_contig_lengths_from_db(conn, host):
    '''
    reads the sequence lengths pf all contigs connected to a particular host from table contigs
    :param conn:
    :param host:
    :return: a python list of sequence lengths
    '''
    c = conn.cursor()
    sql = ''' SELECT sequence_length
                FROM contigs 
               WHERE host=?
            ORDER BY sequence_length DESC'''
    c = conn.cursor()
    rows = c.execute(sql,(host,)).fetchall()
    c.close()
    return rows


def insert_L50_and_contig_num_to_DB(conn, host, number_of_contigs, L50):
    '''
    '''
    c = conn.cursor()
    sql = ''' UPDATE hosts
                 SET number_of_contigs=?, L50=?
               WHERE host=?'''
    c.execute(sql, (number_of_contigs, L50, host))
    c.close()


def add_number_of_contigs_and_L50(database):
    '''
    '''
    conn = connect_to_db(database)
    with conn:
        # Add number of contigs and L50 to the host table
        for host_row in read_hosts_from_db(conn):
            host = host_row[0]
            contig_rows = read_contig_lengths_from_db(conn, host)
            total_seq_length = sum([x[0] for x in contig_rows])
            partial_seq_length = 0
            for i in range(len(contig_rows)):
                partial_seq_length += contig_rows[i][0]
                if partial_seq_length > total_seq_length/2:
                    insert_L50_and_contig_num_to_DB(conn, host, len(contig_rows), i+1)
                    break
    conn.close()


def main(test_flag, database):

    import_gbk_files_to_DB(database, test_flag)
    add_number_of_contigs_and_L50(database)


# Main program ================================================================= 

if __name__ == "__main__":

    print("Commandline input:", " ".join(sys.argv))

    parser = argparse.ArgumentParser()
    parser.add_argument('-db', action="store", dest="database", type=str, default="database.db", help='Name of the sqlite database (default=database.db')
    parser.add_argument('--test', action="store_true", help="creates a test database with only 20 hosts")
    
    args = parser.parse_args()
    database = args.database
    test_flag = args.test

    print("{}: Program started. Database: {}".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), database))
    main(test_flag, database=database)
    print("{}: Program finished".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))