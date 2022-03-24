import sys
import argparse
import sqlite3
from general_functions import connect_to_db, get_base_dir
from datetime import datetime
import os

# Functions ====================================================================


def create_database_tables(database):
    '''
    Creates all necessary tables
    '''
    conn = connect_to_db(database)
    c = conn.cursor()
    c.execute('''
                CREATE TABLE "hosts" (
                	"host"	TEXT NOT NULL UNIQUE,
                	"description"	TEXT NOT NULL,
                	"hosttype"	TEXT NOT NULL,
                	"number_of_contigs"	INTEGER,
                	"L50"	INTEGER,
                	PRIMARY KEY("host")
                )
                '''
    )
    c.execute('''
                CREATE TABLE "contigs" (
	                "host"	TEXT NOT NULL,
	                "contig"	TEXT NOT NULL UNIQUE,
	                "sequence_length"	INTEGER NOT NULL,
	                FOREIGN KEY("host") REFERENCES "hosts"("host"),
	                PRIMARY KEY("contig")
                )
                '''
    )
    c.execute('''
                CREATE TABLE "cds" (
                	"locus_tag"	TEXT NOT NULL UNIQUE,
                	"contig"	TEXT NOT NULL,
                	"product"	TEXT NOT NULL,
                	"translation"	TEXT NOT NULL,
                	"cds_start"	INTEGER NOT NULL,
                	"cds_end"	INTEGER NOT NULL,
                	"core_genome"   TEXT DEFAULT "",
                	FOREIGN KEY("contig") REFERENCES "contigs"("contig"),
                	PRIMARY KEY("locus_tag")
                )
                '''
    )
    c.execute('''
                CREATE TABLE "pfams" (
                	"ID"	INTEGER NOT NULL UNIQUE,
                	"pfamnumber"	TEXT NOT NULL,
                    "locus_tag"    TEXT NOT NULL,
                	"pfamstart"	INTEGER NOT NULL,
                	"pfamend"	INTEGER NOT NULL,
                	"strand"	INTEGER NOT NULL,
	                "contig"	TEXT NOT NULL,
	                FOREIGN KEY("locus_tag") REFERENCES "cds"("locus_tag"),
                    FOREIGN KEY("contig") REFERENCES "contigs"("contig"),
	                PRIMARY KEY("ID" AUTOINCREMENT)
                )
                '''
    )
    c.execute('''
                CREATE TABLE "cluster" (
                	"ID"	INTEGER NOT NULL UNIQUE,
	                "core_genome_indicator"	REAL,
	                "transporter_indicator"	REAL,
	                "number_core_pfams"	INTEGER,
                    "manual_exception" INTEGER,
                	PRIMARY KEY("ID")
                )
                '''
    )
    c.execute('''
                CREATE TABLE "hits" (
                	"hitID"	INTEGER NOT NULL,
                	"clusterID"	INTEGER,
                	"query_contig"	TEXT NOT NULL,
                	"query_pfamID_start"	INTEGER NOT NULL,
                	"query_pfamID_end"	INTEGER NOT NULL,
                	"ref_contig"	TEXT NOT NULL,
                	"ref_pfamID_start"	INTEGER NOT NULL,
                	"ref_pfamID_end"	INTEGER NOT NULL,
                	PRIMARY KEY("hitID" AUTOINCREMENT),
                	FOREIGN KEY("query_pfamID_start") REFERENCES "pfams"("ID"),
                	FOREIGN KEY("ref_pfamID_start") REFERENCES "pfams"("ID"),
                	FOREIGN KEY("ref_contig") REFERENCES "contigs"("contig"),
                	FOREIGN KEY("query_contig") REFERENCES "contigs"("contig"),
                	FOREIGN KEY("ref_pfamID_end") REFERENCES "pfams"("ID"),
                	FOREIGN KEY("query_pfamID_end") REFERENCES "pfams"("ID")
                )
                '''
    )
    c.execute('''
                CREATE TABLE "host_comparisons" (
                	"query_host"	TEXT NOT NULL,
                	"ref_host"	TEXT NOT NULL,
                	FOREIGN KEY("query_host") REFERENCES "hosts"("host"),
                	FOREIGN KEY("ref_host") REFERENCES "hosts"("host"),
                	PRIMARY KEY("query_host","ref_host")
                )
                '''
    )
    
    conn.commit()
    c.close()


def main(database):

    database_path = os.path.join(get_base_dir(), "data", "database", database)
    if os.path.exists(database_path):
        os.remove(database_path)
        print("{}: Already existing database with the same name has been deleted.".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))

    create_database_tables(database)
    print("{}: New database successfully created.".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))


# Main program =================================================================

if __name__ == "__main__":

    print("Commandline input:", " ".join(sys.argv))

    parser = argparse.ArgumentParser()
    parser.add_argument('-db', action="store", dest="database", type=str, default="database.db", help='Name of the sqlite database (default=database.db)')
    
    args = parser.parse_args()
    database = args.database

    print("{}: Program started. Database: {}".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), database))
    main(database=database)
    print("{}: Program finished".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))

