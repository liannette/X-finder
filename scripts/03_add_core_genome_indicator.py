# Consider to include this script into the import_genomes script, as it takes a crazy amount of time to write the core genome information to the database


import sqlite3
import os.path
from subprocess import Popen, PIPE
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import tempfile
import argparse
from general_functions import connect_to_db, timepoint, showtimepoints
import sys
import pandas as pd
import shutil
from datetime import datetime


def get_all_ref_cds(conn):
    """
    """
    conn = connect_to_db(database)
    c = conn.cursor()
    cds_list = c.execute(''' 
                            SELECT locus_tag, translation
                            FROM cds
                        ''') \
                .fetchall()  
    conn.close()  
    return cds_list


def insert_core_genome_in_db(database, core_genome_locus_tags, not_core_genome_locus_tags):
    """
    CDS that are core genome are marked with "-"
    CDS that are not core genome are marked with "+"
    """
    core_genome_locus_tags = [[x] for x in core_genome_locus_tags]
    not_core_genome_locus_tags = [[x] for x in not_core_genome_locus_tags]

    conn = connect_to_db(database)
    with conn:
        c = conn.cursor()
        sql = """ 
                UPDATE cds
                SET core_genome = "-"
                WHERE locus_tag = ? 
              """ 
        c.executemany(sql, core_genome_locus_tags)
        sql = """ 
                UPDATE cds
                SET core_genome = "+"
                WHERE locus_tag = ? 
              """ 
        c.executemany(sql, not_core_genome_locus_tags)
        c.close()
    conn.close()


def main(database):

    dirname = os.path.dirname
    BASE_DIR = dirname(dirname(os.path.abspath(__file__)))
    TEMP_DIR = os.path.join(BASE_DIR, "data", "core genome", "tempdir") # temporary directory. Maybe a dir in memory (tempfile.mkdtemp()) is faster?
    os.mkdir(TEMP_DIR)  

    # Write a fasta file with all reference cds translations
    cds_list = get_all_ref_cds(database)
    print("{}: Fetched all cds from the database".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
    trans_fasta = os.path.join(TEMP_DIR, "translations.fasta")
    with open(trans_fasta, "w") as outfile:
        for locus_tag, translation in cds_list:
            SeqIO.write(SeqRecord(Seq(translation), id=locus_tag, description=""), outfile, "fasta")
    print("{}: Wrote fasta file with all reference cds translations".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))

    # Create diamond databank of the core genome
    diamond_db = os.path.join(TEMP_DIR, "core_genome.dmnd")
    infile = os.path.join(BASE_DIR, "data", "core genome", "Streptomyces_coregenome.fasta") # this is the core genome that was created with orthovenn
    process = Popen(["diamond", "makedb", "--in", infile, "-d", diamond_db], stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()
    print("{}: Created diamond databank".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
    #print("stdout:\n", stdout)
    #print("stderr:\n", stderr)

    # align the cds translation against the core genome
    result_file = os.path.join(TEMP_DIR, "diamond_out.tsv")
    process = Popen(["diamond", "blastp", "-q", trans_fasta , "-d", diamond_db, "-o", result_file, "--fast"], stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()
    print("{}: Diamond finished".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
    #print("stdout:\n", stdout)
    #print("stderr:\n", stderr)

    # Mark the cds that aligned as core genome in the database
    df = pd.read_csv(result_file, sep='\t', header=None)
    core_genome_locus_tags = list(df.iloc[:, 0].unique())
    not_core_genome_locus_tags = [x[0] for x in cds_list if x[0] not in core_genome_locus_tags]
    insert_core_genome_in_db(database, core_genome_locus_tags, not_core_genome_locus_tags)
    print("{}: Inserted core genome information into database".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))

    # Delete temporary directory 
    shutil.rmtree(TEMP_DIR)



if __name__ == "__main__":

    print("Commandline input:", " ".join(sys.argv))

    parser = argparse.ArgumentParser()
    parser.add_argument('-db', action="store", dest="database", type=str, default="database.db", help='Name of the sql database (default=database.db)')
    args = parser.parse_args()
    database = args.database

    print("{}: Program started. Database: {}".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), database))
    main(database=database)
    print("{}: Program finished".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))