import sqlite3
import sys
import argparse
from datetime import datetime
from general_functions import connect_to_db, get_base_dir
import pandas as pd
import math
import os


def get_df_of_hitIDs_and_core_genome(conn):
    """
    """
    sql = """
            WITH distict_cds AS (
            		SELECT DISTINCT hits.hitID, pfams.locus_tag, core_genome
            		  FROM hits
                INNER JOIN pfams ON pfams.pfamID BETWEEN ref_pfamID_start AND ref_pfamID_end
            	INNER JOIN cds   ON pfams.locus_tag = cds.locus_tag
            	)
			  SELECT distict_cds.hitID, 
            	     1.0 * SUM(CASE WHEN distict_cds.core_genome = '-' THEN 1 ELSE 0 END) / COUNT(*) AS core_genome_fraction
			    FROM distict_cds
            GROUP BY distict_cds.hitID
          """
    df = pd.read_sql(sql, conn)  
    df["hitID"] = df["hitID"].map(str)
    df = df.set_index("hitID")
    return df


def get_set_of_transporter_pfams():
    """ 
    Get all pfams associated with transporter activity 
    """
    file_path = os.path.join(get_base_dir(), "data", "transporter pfams", "all_transporter_pfams.txt")
    if os.path.exists(file_path):
        with open(file_path, "r") as infile:
            transporter_pfams = set(infile.read().splitlines())
    else:
        print("{}: File {} not found. Finding all PFAMs that are associated with transporter function and writing them to a file.".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), file_path))
        from get_transporter_pfams import get_transporter_pfams_from_pfam2go
        transporter_pfams = get_transporter_pfams_from_pfam2go()
        with open(file_path, "w") as outfile:
            for pfam in transporter_pfams:
                print(pfam, file=outfile)
    return transporter_pfams


def get_all_clusterIDs(conn):
    '''
    Get a list of all distinct clusterIDs from the hits table
    '''
    sql =   '''   
            SELECT DISTINCT clusterID
            FROM hits 
            ''' 
    c = conn.cursor()
    rows = c.execute(sql).fetchall()
    c.close()
    return [row[0] for row in rows]


def index_hits_clusterID(conn):
    with conn:
        c = conn.cursor()
        c.execute(''' 
                CREATE INDEX hits_clusterID_idx
                ON hits (clusterID)
                ''')
        c.close()


def get_hitIDs_in_cluster(conn, clusterID):
    '''
    '''
    c = conn.cursor()
    sql = """
            SELECT hitID
            FROM hits 
            WHERE clusterID = ?
          """
    rows = c.execute(sql, (clusterID,)).fetchall()
    c.close()
    return [str(row[0]) for row in rows]


def get_core_genome_indicator(conn, hitIDs, core_genome_df):
    """ """
    core_genome_indicator = max(core_genome_df.loc[hitIDs, "core_genome_fraction"])
    return core_genome_indicator


def get_num_core_pfams_and_transporter_indicator(conn, hitIDs_in_cluster, transporter_pfams):
    """
    """
    # get the core pfams of cluster (pfam numbers that are present in each hit)
    core_pfams = None
    for i in range(0, len(hitIDs_in_cluster), 900):  # max number of SQL variables in one query is 999
        hitIDs = hitIDs_in_cluster[i:i+900]
        for host_type in ["ref", "query"]:
            # hit_count specifies in how many sublists the pfam is present
            sql = """
                    WITH distict_pfam AS (
                    SELECT DISTINCT hitID, pfamnumber
                                		  FROM hits
                                    INNER JOIN pfams ON pfams.pfamID BETWEEN {0}_pfamID_start AND {0}_pfamID_end
                    				WHERE hitID IN ({1})
                    				)
                    SELECT pfamnumber, COUNT(*) AS hit_count
                    FROM distict_pfam
                    GROUP BY pfamnumber
                    ORDER BY pfamnumber
                  """.format(host_type, ",".join(hitIDs))
            df = pd.read_sql(sql, conn)
            core_pfams_chunk = set(df[df["hit_count"] == len(hitIDs)]["pfamnumber"]) # take only thos pfams that occured in all sublists
            if core_pfams is None: 
                core_pfams = core_pfams_chunk
            else:
                core_pfams = core_pfams.intersection(core_pfams_chunk)
    # Check how many core_pfams are associated with transporter activity
    if len(core_pfams) == 0:
        fraction_transporter_pfams = 0
    else:
        fraction_transporter_pfams = len(core_pfams.intersection(transporter_pfams))/len(core_pfams)
    return fraction_transporter_pfams, len(core_pfams)


def insert_cluster_to_db(conn, clusterID, core_genome_indicator, transporter_indicator, num_core_pfams):
    """
    """
    with conn:
        c = conn.cursor()
        sql =   """ 
                INSERT OR REPLACE INTO cluster (clusterID, core_genome_indicator, transporter_indicator, number_core_pfams) 
                VALUES (?, ?, ?, ?) 
                """
        c.execute(sql, (clusterID, core_genome_indicator, transporter_indicator, num_core_pfams))
        c.close()


def drop_index_hits_clusterID(conn):
    with conn:
        c = conn.cursor()
        c.execute('''DROP INDEX hits_clusterID_idx''')
        c.close()


def main(database):
    
    conn = connect_to_db(database)

    # Get some general information on core genome and transporter pfams
    core_genome_df = get_df_of_hitIDs_and_core_genome(conn)
    transporter_pfams_set = get_set_of_transporter_pfams()

    # Create index, skip if index already exists. This could be because the program was interruped previously
    try:
        index_hits_clusterID(conn)
    except:
        pass

    clusterIDs = get_all_clusterIDs(conn)
    next_percentage_to_print_status = 0

    for i in range(len(clusterIDs)):
        clusterID = clusterIDs[i]
        hitIDs_in_cluster = get_hitIDs_in_cluster(conn, clusterID)
        core_genome_indicator = get_core_genome_indicator(conn, hitIDs_in_cluster, core_genome_df)
        transporter_indicator, num_core_pfams = get_num_core_pfams_and_transporter_indicator(conn, hitIDs_in_cluster, transporter_pfams_set)
        insert_cluster_to_db(conn, clusterID, core_genome_indicator, transporter_indicator, num_core_pfams)

        # Print progress
        progress_percent = ( (i+1) / len(clusterIDs)) * 100
        if progress_percent > next_percentage_to_print_status:
            print("{}: {}% of cluster finished.".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), math.floor(progress_percent)))
            next_percentage_to_print_status = math.ceil(progress_percent) + 1
    
    drop_index_hits_clusterID(conn)
    conn.close()


# Main program =================================================================

if __name__ == "__main__":

    print("Commandline input:", " ".join(sys.argv))

    parser = argparse.ArgumentParser()
    parser.add_argument('-db', action="store", dest="database", type=str, default="database.db", help='Name of the sql database (default=database.db)')

    args = parser.parse_args()
    database = args.database

    print("{}: Program started. Database: {}".format(
        datetime.now().strftime("%d/%m/%Y %H:%M:%S"), database))
    main(database=database)
    print("{}: Program finished".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))