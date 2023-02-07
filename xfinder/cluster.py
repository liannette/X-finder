import sqlite3
from joblib import Parallel, delayed    # For multithreading
import contextlib
from xfinder.common import print_stdout, print_stderr
import traceback
import sys


def _get_all_sublists(database_path):
    conn = sqlite3.connect(database_path)
    with contextlib.closing(conn.cursor()) as c:
        sql = '''
            SELECT host_type, first_pfamID, last_pfamID 
              FROM sublists
            '''
        rows = c.execute(sql).fetchall()
    conn.close()
    return [tuple(row) for row in rows]



def _find_similar_sublists(database_path, sublist):
    """ 
    Finds similar sublists. 
    A similar sublist is defined by having +/-10% identity of PFAMs.

    Using above definition, sublistA containing 9 PFAMs has only similar
    sublists with the exact same PFAMs. However, for a sublistB
    containing the same PFAMS plus an addional one, sublistA would be a
    similar sublist.
    """
    host_type, first_pos, last_pos = sublist

    conn = sqlite3.connect(database_path)
    with contextlib.closing(conn.cursor()) as c:
        
        # find sublists that form a hit together
        other_host_type = "ref" if host_type  == "query" else "query"
        sql = """
            SELECT "{1}", {1}_first_pfamID, {1}_last_pfamID
              FROM hits
             WHERE {0}_first_pfamID = ?
                   AND {0}_last_pfamID = ?
            """.format(host_type, other_host_type)
        rows = c.execute(sql, (first_pos, last_pos)).fetchall()
        similar_hits = [tuple(row) for row in rows]
        
        # Find sublists that overlap
        sublist_length = last_pos - first_pos + 1  # positions aren't 0 based
        max_pfam_difference = int(sublist_length * 0.1)
        # temporary table number_of_different_pfams contains sublists
        # that have the length of the max_pfam_difference added on both
        # side of sublist. In the next step only the sublists that
        # match the criteria of the max identity difference are chosem
        sql = """
              WITH number_of_different_pfams AS (
                   SELECT host_type, first_pfamID, last_pfamID, 
                   	      ABS(first_pfamID-{0})+ABS(last_pfamID-{1}) 
                              AS pfam_difference
                     FROM sublists
                    WHERE first_pfamID >= {0}-{2}
                          AND last_pfamID <= {1}+{2}
                   )
            SELECT host_type, first_pfamID, last_pfamID
              FROM number_of_different_pfams
             WHERE pfam_difference <= {2}
            """.format(first_pos, last_pos, max_pfam_difference)
        rows = c.execute(sql).fetchall()
        similar_hits += [tuple(row) for row in rows]
        
    conn.close()

    return similar_hits


def cluster_sublists(threads, database_path):

    all_sublists = _get_all_sublists(database_path)

    # Create two dicts, one has the sublist as key and the cluster 
    # number as value, the other has the cluster number as key and a set
    # of all sublists as value
    sublist_cluster_dict = dict.fromkeys(all_sublists, None) 
    cluster_sublists_dict = dict()
    cluster = 0

    for i in range(len(all_sublists)):
        # Check if the sublist is already in a cluster
        if sublist_cluster_dict[all_sublists[i]] is None:

            # search_queue contains sublists that haven't yet been used 
            # to search for similar sublists
            search_queue = set([all_sublists[i]])
            # sublists_in_cluster contains all sublists that have been 
            # already used to search for similar sublists
            sublists_in_cluster = set()
            
            merge_cluster = set()
            
            while len(search_queue) > 0:
                
                # Find similar sublists (matched in a hit or oerlapping)
                results = Parallel(n_jobs=threads) \
                    (delayed(_find_similar_sublists) \
                    (database_path, sublist) for sublist in search_queue)
                similar_sublists = [sublist for r in results for sublist in r]

                # update the sublists in cluster and reset queue
                sublists_in_cluster.update(search_queue)
                search_queue = set()

                for sublist in similar_sublists:
                    # sublist has been clustered before
                    if sublist_cluster_dict[sublist] is not None:
                        # get the cluster to merge
                        merge_cluster.add(sublist_cluster_dict[sublist])
                    # sublist has not been clustered before 
                    elif sublist not in sublists_in_cluster:
                        search_queue.add(sublist)

            # Merge the cluster
            for m_cstr in list(merge_cluster):
                merge_sublists = cluster_sublists_dict.pop(m_cstr)
                sublists_in_cluster.update(merge_sublists)
                
            # Add cluster to dicts
            cluster_sublists_dict[cluster] = sublists_in_cluster
            for sublist in sublists_in_cluster:
                sublist_cluster_dict[sublist] = cluster

            cluster += 1

    # Return a list of clustered hits
    return sorted(list(cluster_sublists_dict.values()), key=len, reverse=True)


def _add_clusterID_to_sublists_table(cluster_id, sublists, conn):
    for sublist in sublists:
        with contextlib.closing(conn.cursor()) as c:
            # Add clusterID to table sublists
            sql =   """
                    UPDATE sublists
                       SET clusterID = ?
                     WHERE first_pfamID = ?
                           AND last_pfamID = ?
                    """
            c.execute(sql, (cluster_id, sublist[1], sublist[2]))


def _get_max_core_genome_fraction(cluster_id, conn):
    with contextlib.closing(conn.cursor()) as c:
        sql = """
            SELECT max(core_genome_fraction)
            FROM sublists
            WHERE clusterID = ?
            """
        c.execute(sql, [str(cluster_id)])
        max_core_genome_fraction = c.fetchall()[0][0]
    return max_core_genome_fraction
    

def _get_antismash_fraction_range(cluster_id, conn):
    with contextlib.closing(conn.cursor()) as c:
        sql = """
            SELECT min(antismash_fraction), max(antismash_fraction)
            FROM sublists
            WHERE clusterID = ?
            """
        c.execute(sql, [str(cluster_id)])
        min_as_fraction, max_as_fraction = c.fetchall()[0]
    return min_as_fraction, max_as_fraction


def _get_transporter_indicator(cluster_id, sublists, conn):
    """
    """
    with contextlib.closing(conn.cursor()) as c:
        sql = """
            WITH pfam_count AS (
                SELECT pfam_num, transporter, 
                       COUNT(DISTINCT sublists.ROWID) as cnt
                  FROM sublists
                       INNER JOIN pfams 
                       ON pfamID BETWEEN first_pfamID AND last_pfamID
                 WHERE clusterID = ?
			     GROUP by pfam_num
           		)
            SELECT AVG(transporter), COUNT(*)
            FROM pfam_count
            WHERE cnt = ?
            """
        c.execute(sql, [cluster_id, len(sublists)])
        transporter_ind, number_core_pfams = c.fetchall()[0]
        if transporter_ind is None:
            transporter_ind = 0
        
    return round(transporter_ind, 3), number_core_pfams


def _add_cluster_information(cluster_id, max_core_genome_fraction, 
                             antismash_fraction_range, transporter_ind, 
                             number_core_pfams, conn):

    with contextlib.closing(conn.cursor()) as c:
        sql = '''
            INSERT INTO cluster (clusterID, 
                                 max_core_genome_fraction, 
                                 min_antismash_fraction, 
                                 max_antismash_fraction,
                                 transporter_indicator, 
                                 number_core_pfams)
                 VALUES (?, ?, ?, ?, ?, ?)
            '''
        c.execute(sql, (cluster_id, 
                        max_core_genome_fraction,
                        antismash_fraction_range[0], 
                        antismash_fraction_range[1], 
                        transporter_ind, 
                        number_core_pfams
                        ))


def add_cluster_to_db(clustered_sublists, database_path):
    conn = sqlite3.connect(database_path)
    # # Auto-roll back if a sql error occurs
    # with conn:
    for i in range(len(clustered_sublists)):
        cluster_id = i+1
        sublists = list(clustered_sublists[i])
        _add_clusterID_to_sublists_table(cluster_id, sublists, conn)
        conn.commit()
        max_core_genome_fraction = _get_max_core_genome_fraction(cluster_id, 
                                                                 conn)
        antismash_fraction_range = _get_antismash_fraction_range(cluster_id, 
                                                                 conn)
        transporter_ind, number_core_pfams = _get_transporter_indicator(
            cluster_id, sublists, conn)
        _add_cluster_information(
            cluster_id, 
            max_core_genome_fraction, 
            antismash_fraction_range, 
            transporter_ind, 
            number_core_pfams, 
            conn)
    conn.close()


def cluster_hits(threads, database_path, out_dir):
    try:
        print_stdout("Grouping similar hits into clusters", out_dir)
        clustered_sublists = cluster_sublists(threads, database_path)
        add_cluster_to_db(clustered_sublists, database_path)
        print_stdout(f"Hits were successfully grouped into "
                     f"{len(clustered_sublists)} clusters.", out_dir)
    except:
        print_stderr(traceback.format_exc(), out_dir)
        sys.exit(1)