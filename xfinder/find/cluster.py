import sqlite3
from joblib import Parallel, delayed    # For multithreading
import contextlib



def _get_all_hit_ids(conn):
    with contextlib.closing(conn.cursor()) as c:
        rows = c.execute('''SELECT hitID FROM hits''').fetchall()
    return [row[0] for row in rows]


def _hit_to_sublists(hit):
    hitID = hit[0]
    query_sublist = ("query", *hit[1:3])
    ref_sublist = ("ref", *hit[3:5])
    return hitID, query_sublist, ref_sublist


def _get_sublists_from_db(conn, hit_id):
    """
    Gets the first and the last pfamID for both query and reference
    sublist of a hit
    """
    with contextlib.closing(conn.cursor()) as c:
        sql =   """
                SELECT hitID, 
                       query_first_pfamID, query_last_pfamID, 
                       ref_first_pfamID, ref_last_pfamID
                FROM hits
                WHERE hitID = ?
                """
        row = c.execute(sql, (hit_id,)).fetchall()[0]
    _, query_sublist, ref_sublist = _hit_to_sublists(row)
    return query_sublist, ref_sublist


def _find_similar_hits(database_path, sublist):
    """ 
    Finds similar sublists. 
    A similar sublist is defined by having +/-10% identity of PFAMs.

    Using above definition, sublistA containing 9 PFAMs has only similar
    sublists with the exact same PFAMs. However, for a sublistB
    containing the same PFAMS plus an addional one, sublistA would be a
    similar sublist.

    If a sublist already has a clusterID, if will not be put into the
    queue, but the clusters will be merged later.
    """
    host_type, first_pos, last_pos = sublist
    sublist_length = last_pos - first_pos + 1 # +1 as positions aren't 0 based
    max_pfam_difference = int(sublist_length * 0.1)

    conn = sqlite3.connect(database_path)
    with contextlib.closing(conn.cursor()) as c:
        # temporary table number_of_different_pfams contains only hits
        # that have the length of the max_pfam_difference added on both
        # side ofthe ref/query pfam sublist. In the next step only the
        # hits are taken that matches the criteria of the max identity
        # difference
        sql = """
              WITH number_of_different_pfams AS (
                   SELECT hitID,
                          query_first_pfamID, query_last_pfamID, 
                          ref_first_pfamID, ref_last_pfamID,
                   	      ABS({0}_first_pfamID-{1})+ABS({0}_last_pfamID-{2}) 
                              AS pfam_difference
                     FROM hits
                    WHERE {0}_first_pfamID >= {1}-{3}
                          AND {0}_last_pfamID <= {2}+{3}
                   )
            SELECT hitID,
                   query_first_pfamID, query_last_pfamID, 
                   ref_first_pfamID, ref_last_pfamID
              FROM number_of_different_pfams
             WHERE pfam_difference <= {3}
            """.format(host_type, first_pos, last_pos, max_pfam_difference)
        rows = c.execute(sql).fetchall()
    conn.close()

    return rows


def cluster_hits(conn, threads, database_path):

    conn = sqlite3.connect(database_path)

    all_hit_ids = _get_all_hit_ids(conn)
    hit_to_cluster = dict.fromkeys(all_hit_ids, None) # cluster is value
    cluster_to_hits = dict() # cluster is key, set of hitIDs is value
    cluster = 0

    for i in range(len(all_hit_ids)):
        hit_id = all_hit_ids[i]
        if hit_to_cluster[hit_id] is None:
            
            # sublists_queue contains sublists that are in the cluster, 
            # but haven't yet been used to search for similar sublists
            sublists_queue = set(_get_sublists_from_db(conn, hit_id))
            # sublists_in_cluster contains all sublists that have been 
            # already used to search for similar sublists
            sublists_in_cluster = set()
            hit_ids_in_cluster = set()

            while len(sublists_queue) > 0:
                
                # Find hits with similar sublists
                results = Parallel(n_jobs=threads)(delayed(_find_similar_hits) \
                    (database_path, sublist) for sublist in sublists_queue)
                similar_hits = [hit for r in results for hit in r]

                sublists_in_cluster.update(sublists_queue)
                sublists_queue = set() # reset queue

                for hit in similar_hits:
                    hit_id, query_sublist, ref_sublist = _hit_to_sublists(hit)
                    if hit_to_cluster[hit_id] is None:
                        # Hit has not been put into a cluster before
                        hit_ids_in_cluster.add(hit_id)
                        for sublist in (query_sublist, ref_sublist):
                            if sublist not in sublists_in_cluster:
                                sublists_queue.add(sublist)
                    else:
                        # Hit is has been clustered before
                        merge_cluster = hit_to_cluster[hit_id]
                        merge_hit_ids = cluster_to_hits.pop(merge_cluster)
                        hit_ids_in_cluster.update(merge_hit_ids)

            # Add cluster information to dicts
            cluster_to_hits[cluster] = hit_ids_in_cluster
            for hit_id in hit_ids_in_cluster:
                hit_to_cluster[hit_id] = cluster

            cluster += 1
        
    conn.close()

    # Return a list of clustered hits
    return sorted(list(cluster_to_hits.values()), key=len, reverse=True)


def add_cluster_information_to_db(clustered_hits, conn):
    """
    """
    with conn:
        with contextlib.closing(conn.cursor()) as c:

            for i in range(len(clustered_hits)):
                cluster_id = i+1
                hit_ids_in_cluster = list(clustered_hits[i])

                # max number of SQL variables (?) in one query is 999
                for j in range(0, len(hit_ids_in_cluster), 990):
                    hit_ids = hit_ids_in_cluster[j:j+900]
                    sql =   """
                            UPDATE hits
                            SET clusterID = ?
                            WHERE hitID in ({})
                            """.format(','.join(['?']*len(hit_ids)))
                    c.execute(sql, (cluster_id, *hit_ids))


#------------------------------------------------------------------------------

import pandas as pd

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

def get_core_genome_indicator(hit_ids, core_genome_df):
    """ """
    core_genome_indicator = max(core_genome_df.loc[hit_ids, "core_genome_fraction"])
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