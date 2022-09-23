import sqlite3
import sys
from general_functions import connect_to_db, get_base_dir
import argparse
import pandas as pd
import os
from datetime import datetime
from joblib import Parallel, delayed    # For multithreading
import math


def index_sublist_positions(conn):
    with conn:
        c = conn.cursor()
        c.execute(''' 
                CREATE INDEX query_sublist_pos
                ON hits (query_pfamID_start, query_pfamID_end) 
                ''')
        c.execute(''' 
                CREATE INDEX ref_sublist_pos
                ON hits (ref_pfamID_start, ref_pfamID_end) 
                ''')
        c.close()


def get_highest_hitID(conn):
    with conn:
        c = conn.cursor()
        max_hitID = c.execute('''SELECT MAX(hitID) FROM hits''').fetchall()[0][0]
        c.close()
    return max_hitID


def get_sublists(conn, hitID):
    """
    Gets the position of the first and the last PFAM number 
    in regards to the pfam genome sequence table.
    """
    c = conn.cursor()
    sql =   """
            SELECT query_pfamID_start, query_pfamID_end, ref_pfamID_start, ref_pfamID_end
            FROM hits
            WHERE hitID = ?
            """
    rows = c.execute(sql, (hitID,)).fetchall()
    sublist_positions = hit_row_to_sublist_positions(rows[0])
    return sublist_positions


def hit_row_to_sublist_positions(row):
    return ("query", *row[0:2]) , ("ref", *row[2:4])


def find_similar_hits(database, sublist_pos):
    """ 
    We also need to check if a hit is already in a cluster. Imagine following case: 
    A previous clustered hit contained a sublist with start 100 and end 109. 
    If a similar sublist is defined by +/-10% identity it would clusters only 
    with identical sublists. But a if later on there is a hit with start 100 and end 
    110, it would cluster with the previous one. Therefore the cluster must be combined.
    If a sublist already has a clusterID, if will not put into the queue, but the clusterID
    will be used later on to add all sublists belonging to that clusterID into the new
    cluster.
    """
    conn = connect_to_db(database)
    c = conn.cursor()
    hosttype, start_pos, end_pos = sublist_pos
    sql = """
        WITH number_diff_pfams AS (
            SELECT hitID, clusterID, query_pfamID_start, query_pfamID_end, ref_pfamID_start, ref_pfamID_end,
            	   ABS({0}_pfamID_start-{1})+ABS({0}_pfamID_end-{2}) AS diff
            FROM hits
            WHERE {0}_pfamID_start >= {1}-{3}
              AND {0}_pfamID_end <= {2}+{3}
            )
        SELECT query_pfamID_start, query_pfamID_end, ref_pfamID_start, ref_pfamID_end, hitID, clusterID
        FROM number_diff_pfams
        WHERE diff <= {3}
        """.format(hosttype, start_pos, end_pos, (end_pos-start_pos+1)/10)
    rows = c.execute(sql).fetchall()
    conn.close()

    hitIDs = set()
    sublists_queue = set()
    clusterIDs = set()
    for row in rows:
        if row[5] is None:
            hitIDs.add(row[4])
            sublists_queue.update(hit_row_to_sublist_positions(row))
        else:
            clusterIDs.add(row[5])

    return hitIDs, sublists_queue, clusterIDs


def update_hits_with_clusterID(conn, clusterID, hitIDs_in_cluster, old_clusterIDs):
    """
    Updates the table hits by adding the clusterID
    """
    with conn:
        c = conn.cursor()
        hitIDs_in_cluster = sorted(hitIDs_in_cluster)
        for i in range(0, len(hitIDs_in_cluster), 900):  # max number of SQL variables in one query is 999
            hitIDs = hitIDs_in_cluster[i:i+900]
            sql =   """
                    UPDATE hits
                    SET clusterID = ?
                    WHERE hitID in ({})
                    """.format(','.join(['?']*len(hitIDs)))
            c.execute(sql, (clusterID, *hitIDs))
        # Below a s eperate query for the clusterID, because hitID has a index but clusterID not. 
        # If they would be combined into one query, the "WHERE hitID in" would not use index.
        for i in range(0, len(old_clusterIDs), 900):
            clusterIDs = old_clusterIDs[i:i+900]
            sql =   """
                    UPDATE hits
                    SET clusterID = ?
                    WHERE clusterID in ({})
                    """.format(','.join(['?']*len(clusterIDs)))
            c.execute(sql, (clusterID, *clusterIDs))
        c.close()


def drop_index_sublist_positions(conn):
    with conn:
        c = conn.cursor()
        c.execute(''' DROP INDEX query_sublist_pos ''')
        c.execute(''' DROP INDEX ref_sublist_pos ''')
        c.close()


def cluster_hits(threads, conn):
    
    # Create index, skip if index already exists. This could be because the program was interruped previously
    try:
        index_sublist_positions(conn)
    except:
        pass

    cluster_cnt = 1
    already_clustered_hitIDs = set()
    next_percentage_to_print_status = 0
    num_of_hits = get_highest_hitID(conn)

    for hitID in range(1, num_of_hits+1):

        if hitID not in already_clustered_hitIDs:

            sublists_queue = set(get_sublists(conn, hitID)) # Contains sublist positions that has not been used yet to query the database
            sublists_in_cluster = set().union(sublists_queue)
            hitIDs_in_cluster = set([hitID])
            old_clusterIDs = set()  # If similar sublists have been put into a cluster previously, the clusterID will we stored here

            # collect all hitIDs
            while len(sublists_queue) > 0:
                
                results = Parallel(n_jobs=threads)(delayed(find_similar_hits) \
                    (database, sublist_pos) for sublist_pos in sublists_queue)
                hitIDs, sublists, clusterIDs = zip(*results)
                hitIDs_in_cluster = hitIDs_in_cluster.union(*hitIDs)
                sublists_queue = set().union(*sublists).difference(sublists_in_cluster)
                sublists_in_cluster = sublists_in_cluster.union(sublists_queue)
                old_clusterIDs = old_clusterIDs.union(*clusterIDs)
        
            # Get clusterID
            old_clusterIDs = sorted(old_clusterIDs)
            if len(old_clusterIDs) == 0:
                clusterID = cluster_cnt
                cluster_cnt += 1
            else:
                clusterID = old_clusterIDs.pop(0)

            update_hits_with_clusterID(conn, clusterID, hitIDs_in_cluster, old_clusterIDs)
            already_clustered_hitIDs = already_clustered_hitIDs.union(hitIDs_in_cluster)

            # Print progress
            percent_clustered_hits = (len(already_clustered_hitIDs) / num_of_hits) * 100
            if percent_clustered_hits > next_percentage_to_print_status:
                print("{}: {}% of hits clustered.".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), math.floor(percent_clustered_hits)))
                next_percentage_to_print_status = math.ceil(percent_clustered_hits)
        
    drop_index_sublist_positions(conn)


def main(threads, database):
    
    conn = connect_to_db(database)
    cluster_hits(threads, conn)
    conn.close()


# Main program =================================================================

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