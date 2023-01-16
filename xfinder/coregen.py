import sqlite3
from subprocess import Popen, PIPE
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import pandas as pd
import contextlib
from xfinder.common import print_stdout, print_stderr
import sys
import traceback
import tempfile
import os


def _get_cds(database_path):
    """
    database_path must be absolute path
    """
    conn = sqlite3.connect(database_path)
    with contextlib.closing(conn.cursor()) as c:
        sql = ''' SELECT locus_tag, translation FROM cds '''
        cds_list = c.execute(sql).fetchall()  
    conn.close()  
    return cds_list


def cds_translations_to_fasta(database_path, cds_trans_fasta):
    """ Write a fasta file with all cds translations """
    cds_list = _get_cds(database_path)
    with open(cds_trans_fasta, "w") as outfile:
        for locus_tag, translation in cds_list:
            seq_record = SeqRecord(Seq(translation), id=locus_tag)
            SeqIO.write(seq_record, outfile, "fasta")
    return cds_list
    

def create_diamond_database(diamond_db, core_genome_path, out_dir, threads):
    """ Create diamond databank of the core genome """
    process = Popen(["diamond", "makedb", "--in", core_genome_path, "-d", 
                      diamond_db, "--threads", str(threads)], 
                    stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate() 

    # No error occured
    if process.returncode == 0:
        # Both stdout & stderr are send to stdout. Diamond
        # sends status messages to stderr, even if there is no error
        for msg in (stdout, stderr):
            if len(msg) > 0:
                print_stdout(msg.decode(), out_dir)
    # Error occured      
    else:
        # stderr is send to stderr
        if len(stdout) > 0:
            print_stdout(stdout.decode(), out_dir)
        if len(stderr) > 0:
            print_stderr(stderr.decode(), out_dir)
        raise RuntimeError("Error occured while creating the diamond " 
                           "database.")


def run_diamond(trans_fasta, diamond_db, result_file, out_dir, threads):
    """ align the cds translations against the core genome """
    # Add the number of threads!
    process = Popen(["diamond", "blastp", "-q", trans_fasta , "-d", diamond_db,
                     "-o", result_file, "--threads", str(threads), "--fast"], 
                    stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()

    if process.returncode == 0:
        # Both stdout and stderr are send to stdout. Diamond
        # sends status messages to stderr, even if there is no error
        for msg in (stdout, stderr):
            if len(msg) > 0:
                print_stdout(msg.decode(), out_dir)
    else:
        # Error occured, therefore stderr is also send to stderr
        if len(stdout) > 0:
            print_stdout(stdout.decode(), out_dir)
        if len(stderr) > 0:
            print_stderr(stderr.decode(), out_dir)
        raise RuntimeError("Error occured while running diamond.")


def _core_genome_information_to_db(database_path, core_genome_locus_tags, 
                                   not_core_genome_locus_tags):
    """
    CDS that have aligned to core genome are marked with 1,
    CDS that have not aligned to core genome are marked with 0
    """
    conn = sqlite3.connect(database_path)
    with conn:
        with contextlib.closing(conn.cursor()) as c:
            sql = """ 
                    UPDATE cds
                       SET core_genome = 1
                     WHERE locus_tag = ? 
                  """ 
            c.executemany(sql, [[x] for x in core_genome_locus_tags])
            sql = """ 
                    UPDATE cds
                       SET core_genome = 0
                     WHERE locus_tag = ? 
                  """ 
            c.executemany(sql, [[x] for x in not_core_genome_locus_tags])
    conn.close()


def add_core_genome_information(database_path, result_file, cds_list):
    """
    """
    df = pd.read_csv(result_file, sep='\t', header=None)
    core_genome_locus_tags = set(df.iloc[:, 0])
    not_core_genome_locus_tags = set(
        [x[0] for x in cds_list if x[0] not in core_genome_locus_tags])
    _core_genome_information_to_db(database_path, core_genome_locus_tags, 
                                   not_core_genome_locus_tags)
    return len(core_genome_locus_tags)


def add_coregenome_info(database_path, core_genome_path, out_dir, threads):

    try:
        print_stdout("Adding core genome information. Core genome file: "
                    f"{core_genome_path}", out_dir)

        with tempfile.TemporaryDirectory() as temp_dir:
            
            trans_fasta = os.path.join(temp_dir, "translations.fasta")
            diamond_db = os.path.join(temp_dir, "core_genome.dmnd")
            result_file = os.path.join(temp_dir, "diamond_out.tsv")

            # Write a fasta file with all cds translations
            cds_list = cds_translations_to_fasta(database_path, trans_fasta)
            # Create diamond databank of the core genome
            create_diamond_database(diamond_db, core_genome_path, out_dir, 
                                    threads)
            # Run diamond
            run_diamond(trans_fasta, diamond_db, result_file, out_dir, threads)

            # Add the core genome information to the database
            num_core_cds = add_core_genome_information(
                database_path, result_file, cds_list)
        
        print_stdout("Inserted core genome information into database. "
                    f"{num_core_cds}/{len(cds_list)} CDS translations aligned " 
                    "to the core genome with >90% identity", out_dir)
    except:
        print_stderr(traceback.format_exc(), out_dir)
        sys.exit(1)