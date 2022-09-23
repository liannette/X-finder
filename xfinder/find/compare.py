# Original idea and code from Xinglin Jiang <xinji@biosustain.dtu.dk>
# clearned up with the help of Simon Shaw <sisha@biosustain.dtu.dk>
# rewritten and more features added by Annette Lien <a.lien@posteo.de>

# - one fragment from either of the two sequence can have more than one
# matches on the other, and all of them will be reported.
# - we assume all genomes are in linear form. this may overlook the
# genes near the breaking point of circular genomes.it should not be a
# problem if circular genomes are always break at the ori area. because
# the genes in ori area are already well studied.
# - when pfam have an overlap of at least 50%, only the pfam with the 
# smalles pfam number is taken

import sqlite3
import contextlib
import os.path
import numpy as np
from xfinder.find.common import get_database_dir, print_stderr
from joblib import Parallel, delayed
import itertools
import traceback


class Pfam:

    def __init__(self, pfam_id, pfam_num, cds_locus_tag, start, end, strand, 
                 antismash_core):
        self.pfam_id = pfam_id
        self.pfam_num = pfam_num
        self.cds_locus_tag = cds_locus_tag
        self.start = start  # 0-based
        self.end = end      # 0-based
        self.strand = strand
        self.antismash_core = antismash_core


class SeqRecord:

    def __init__(self, acc, g_pfam_list):
        self.seq_acc = acc
        self.g_pfam_list = g_pfam_list # genome pfam sequence

    @classmethod
    def _from_db(cls, seq_acc, conn):
        ''' Returns the genome pfam sequence of a sequence record '''
        with contextlib.closing(conn.cursor()) as c:
            sql = ''' SELECT pfamID, pfam_num, locus_tag, pfam_start, pfam_end,
                             strand, antismash_core
                        FROM pfams 
                       WHERE seq_acc=? '''
            rows = c.execute(sql,(seq_acc,)).fetchall()
        return cls(acc=seq_acc, 
                   g_pfam_list=[Pfam(*row) for row in rows])
                   
    def _get_pfam_num_seq(self):
        return [pfam.pfam_num for pfam in self.g_pfam_list]

class HitParameters:

    def __init__(self, seed_size, gap_threshold, size_threshold, 
                 dna_length_threshold):
        self.seed_size = seed_size
        self.gap_threshold = gap_threshold
        self.size_threshold = size_threshold
        self.dna_length_threshold = dna_length_threshold


class Hit:

    def __init__(self, query_start, query_end, ref_start, ref_end, 
                 hit_direction):
        if query_start >= query_end:
            raise ValueError(f'query_start must be < query_end')
        if ref_start >= ref_end:
            raise ValueError(f'ref_start must be < ref_end')
        self.query_start = query_start
        self.query_end = query_end
        self.ref_start = ref_start
        self.ref_end = ref_end
        self.hit_direction = hit_direction  # direct==1, inverted==-1

    def get_query_length(self):
        return self.query_end - self.query_start

    def get_ref_length(self):
        return self.ref_end - self.ref_start

    def has_same_query_location(self, other):
        return (self.query_start == other.query_start
                and self.query_end == other.query_end)

    def is_inside_query(self, other):
        return (self.query_start >= other.query_start 
                and self.query_end <= other.query_end)

    def query_overlaps(self, other):
        return (self.query_start < other.query_end 
                and self.query_end > other.query_start)

    def ref_overlaps(self, other):
        return (self.ref_start < other.ref_end 
                and self.ref_end > other.ref_start)

    def both_regions_overlap(self, other):
        return (self.query_overlaps(other) 
                and self.ref_overlaps(other))

    def gap_calculation(self, other):
        '''
        calculates the gap between two hits. a gap value below 0 
        indicates that the two sublist overlap with eachother.
        gaps within the two sublists (from previous combination of hits) 
        were not considered.
        Returns: gap_in_query, gap_in_reference
        '''
        # here we calculate the distance/gap between the two sublists.
        # gaps within the two sublists were not considered
        gap_in_query = max(other.query_start-self.query_end, 
                           self.query_start-other.query_end)
        gap_in_reference = max(other.ref_start-self.ref_end, 
                               self.ref_start-other.ref_end)
        return (gap_in_query, gap_in_reference)

    def combine_hits(self, other):
        '''
        combines two hits into one.
        the direction will be the sum of the directions of all the
        matched pfams don't try to add gap value in to the sublist
        tuple. when two sublist overlap with eachother, there is no way
        to calculate the gap value in the combined_sublist, as that
        require the position of the gaps.
        '''
        combined_start_in_query = min(other.query_start, self.query_start)
        combined_end_in_query = max(other.query_end, self.query_end)
        combined_start_in_reference = min(other.ref_start, self.ref_start)
        combined_end_in_reference = max(other.ref_end, self.ref_end)
        combined_direction = other.hit_direction + self.hit_direction
        combined = Hit(combined_start_in_query, combined_end_in_query,
                       combined_start_in_reference, combined_end_in_reference, 
                       combined_direction)
        return combined

    def __repr__(self):
        rep = ("Hit(" 
               + ", ".join([self.query_start, self.query_end, self.ref_start,
                            self.ref_end, self.hit_direction]) 
               + ")")
        return rep


def count_comparisons(conn):
    '''
    Counts how many entries there are in the table host_comparisons
    '''
    with contextlib.closing(conn.cursor()) as c: # auto-closes cursor
        sql = ''' SELECT COUNT(*) FROM host_comparisons '''
        c.execute(sql,)
        num_comparisons = c.execute(sql).fetchall()[0][0]
    return num_comparisons


def _get_host_ids(conn, host_type, max_l50):
    '''
    Gets the hostIDs for either query or reference hosts
    '''
    with contextlib.closing(conn.cursor()) as c:
        sql = ''' SELECT hostID
                    FROM hosts 
                WHERE host_type = ? AND L50 <= ?'''
        rows = c.execute(sql, (host_type, max_l50)).fetchall()
        hostIDs = [row[0] for row in rows]
    return hostIDs


def get_host_comparison_pairs(conn, max_l50):
    """
    Returns a tuple containing:
    - a list of all host comparisons pairs (query_hostID, ref_hostID) 
    - the number of query hosts 
    - the number of reference hosts.
    """
    ref_host_ids = _get_host_ids(conn, "ref", max_l50)
    query_host_ids = _get_host_ids(conn, "query", max_l50)
    host_pairs = list(itertools.product(query_host_ids, ref_host_ids))
    return host_pairs, len(query_host_ids), len(ref_host_ids)


def _temp_database_path(database_path, thread):
    """ Gets the database path for a thread """
    temp_db =  database_path.split('.')[-2] + '_temp' + str(thread) + '.db'
    return os.path.join(get_database_dir(), temp_db)


def _host_pairs_into_chunks(host_pairs, threads):
    """
    Splits a list of host comparison pairs into equally sized chunks, 
    one for each thread. Return a enumerated list of chunks, so that 
    the chunks can be assigned to a specific thread.
    """
    return list(enumerate(np.array_split(host_pairs, threads)))


def _create_temporary_databases(database_path, threads):
    '''
    Create temporary databases to write the hits with multithreading. 
    Sqlite does not allow concurrent database writes to the same database
    The temporary databases will be merged later into the main one 
    afterwards.
    '''
    for thread in range(threads):
        temp_database_path = _temp_database_path(database_path, thread)
        conn = sqlite3.connect(temp_database_path)
        with conn:
            with contextlib.closing(conn.cursor()) as c:
                c.execute('''
                    CREATE TABLE "hits" (
                	    "query_seq_acc"	        TEXT NOT NULL,
                	    "query_first_pfamID"    INTEGER NOT NULL,
                	    "query_last_pfamID"	    INTEGER NOT NULL,
                	    "ref_seq_acc"	        TEXT NOT NULL,
                	    "ref_first_pfamID"	    INTEGER NOT NULL,
                	    "ref_last_pfamID"	    INTEGER NOT NULL,
                        PRIMARY KEY("query_first_pfamID", "ref_first_pfamID")
                        )
                        '''
                )
                c.execute('''
                        CREATE TABLE "host_comparisons" (
                            "query_hostID"	INTEGER NOT NULL,
                            "ref_hostID"    INTEGER NOT NULL,
                            PRIMARY KEY("query_hostID","ref_hostID")
                        )
                        '''
                )            
        conn.close()


def _host_comparison_exists(conn, query_hostID, ref_hostID):
    '''
    Checks if hits have already been calculated for the two hits
    '''
    with contextlib.closing(conn.cursor()) as c:
        sql = ''' SELECT COUNT(*) 
                    FROM host_comparisons 
                   WHERE query_hostID=? 
                     AND ref_hostID=?
              '''
        result = c.execute(sql,(query_hostID, ref_hostID)).fetchall()[0][0]
    return result == 1


def _get_seq_records(conn, hostID):
    '''
    Returns all sequence accessions belonging to a specific hostID
    '''
    with contextlib.closing(conn.cursor()) as c:
        sql = ''' SELECT seq_acc 
                    FROM seq_records 
                   WHERE hostID=?'''
        rows = c.execute(sql,(hostID,)).fetchall()
        seq_records = [SeqRecord._from_db(row[0], conn) for row in rows]
    return seq_records


def _initialize_2d_array(query_seq, ref_seq):
    """ 
    2D array with query sequence on the X axis and reference on the
    Y axis. "+2", because we need additional rows/columns with
    value 0 at both ends of each sequence
    """
    return np.zeros((2+len(query_seq), 2+len(ref_seq)), dtype=int)


def _direct_match_end(x, y, query_seq, ref_seq, score, length_threshold):
    """ 
    Returns True if (bottom or right of the array is reached or
    continued match stops) and match length >= threshold.
    Otherwise returns False
    """
    if ((x == len(query_seq) - 1
            or y == len(ref_seq) - 1
            or query_seq[x + 1] != ref_seq[y + 1]) 
        and score[x + 1][y + 1] >= length_threshold):
        return True
    else:
        return False


def _inverted_match_end(x, y, query_seq, ref_seq, score_i, length_threshold):
    """
    Returns True if (bottom or left of the array is reached or continued
    match stops) and match length >= threshold.
    """
    if ((x == len(query_seq) - 1
            or y == 0
            or query_seq[x + 1] != ref_seq[y - 1]) 
        and score_i[x + 1][y + 1] >= length_threshold):
        return True
    else:
        return False


def _find_common_sublists(query_seq, ref_seq, length_threshold):
    '''
    purpose: compare a query gpfamnumbersequence with a reference 
    gpfamnumbersequence to find out all common pfamnumber sublists 
    longer than the threshold (2 by default) in both directions.
    the common sublists are all exact matches. insertion, deletion or
    invertion is not tolerated in this step. these represent the exactly
    matched pfam domain sequences from both genomes. pfamnumber repeats 
    (for example from pks and nrps genes) generate overlapped 
    common_sublists. they will be combined in the next step.
    input: a list of pfamnumbers representing the query genome and a
    list represent the reference genomes. a length_threshold of the
    match, it is 2 by default.
    output: a list of hits. 
    a match with length of 2 has a direction of 1 or -1
    '''

    # One arrays for detection of direct (score) and one for inverted
    # matches. the arrays have two extra rows/columns at each side
    score = _initialize_2d_array(query_seq, ref_seq) 
    score_i = _initialize_2d_array(query_seq, ref_seq)
    hits = []
    # classic Dynamic programming for finding common substrings
    for x in range(0, len(query_seq)):
        for y in range(0, len(ref_seq)):

            if query_seq[x] == ref_seq[y]:
                score[x + 1][y + 1] = score[x][y] + 1
                score_i[x + 1][y + 1] = score_i[x][y + 2] + 1
                if _direct_match_end(x, y, query_seq, ref_seq, score, 
                                          length_threshold):
                    hit = Hit(query_start = x + 1 - score[x + 1][y + 1], 
                              query_end = x + 1, 
                              ref_start = y + 1 - score[x + 1][y + 1], 
                              ref_end = y + 1, 
                              hit_direction = score[x + 1][y + 1] - 1)
                    hits.append(hit)
                if _inverted_match_end(x, y, query_seq, ref_seq, score_i, 
                                       length_threshold):
                    hit = Hit(query_start = x + 1 - score_i[x + 1][y + 1],
                              query_end = x + 1,
                              ref_start = y,
                              ref_end = y + score_i[x + 1][y + 1], 
                              hit_direction = -score_i[x + 1][y + 1] + 1)
                    hits.append(hit)
    return hits


def _group_neighboring_hits(hits, gap_threshold):
    '''
    https://en.wikipedia.org/wiki/Approximate_string_matching
    purpose: 1) to combine the sublists that are close to (including
    overlap with) eachother on BOTH the qurery and reference genomes.
    this is similar to but more complex than the classic "merge
    overlapping intervals found in a list" problem
    2)to group the sublists that are close to (including overlap with)
    eachother on both the qurery and reference genomes. one such group
    will represent a gene cluster with some degree of modification,
    such as insertion deletion or invertion. a group is a list of
    sublists. Overlapped sublists generated from pfamnumber repeats
    (for example from pks and nrps genes which have tandom repeats of
    pfam domains) will be kept in the group. The one pfam to one pfam
    corresponding relationship in a group maybe used in futher for
    generating figures of gene cluster alignment. However a group can be
    too complex to be interpreted by human eyes.
    3)gap_threshold is the biggest gap allowed.
    '''
    all_combined_hits = []
    while len(hits) > 0:
        
        combined_hit = hits.pop() # starting with last element in list
        i = len(hits) - 1
        while i >= 0:
            gap_in_query, gap_in_reference \
                    = combined_hit.gap_calculation(hits[i])
            if gap_in_query > gap_threshold:
                i = -1
                # combined_hit is finished. query_gap will only grow for
                # previous elements, as hits are ordered by increasing
                # query_end values
            elif gap_in_reference > gap_threshold:
                i -= 1
                # Compare with next hit 
            else:
                # gap is small enough in both query and reference
                combined_hit = combined_hit.combine_hits(hits.pop(i))
                i = len(hits) - 1
                # Combine hits and start comparing with all hits again
        
        all_combined_hits.append(combined_hit)

    all_combined_hits.reverse()
    return all_combined_hits


def _completely_group_neighboring_hits(hits, gap_threshold):
    """
    Sometimes an already grouped hit can group with another already
    grouped hit.
    """
    while True:
        list_length = len(hits)
        hits = _group_neighboring_hits(hits, gap_threshold)
        if len(hits) == list_length:
            # no grouping happened in this iteration
            break
    return hits


def _filter_by_size(hits, size_threshold):
    """ 
    Only keeps those hits, where the DNA sequence length of both, query 
    and reference sublist, is long enough
    """
    filtered_hits = []
    for hit in hits:
        if (hit.get_query_length() >= size_threshold 
            and hit.get_ref_length() >= size_threshold):
            # combined_sublist is not too short on both of the genomes
            filtered_hits.append(hit)
    return filtered_hits


def _belong_similarity(x,y):
    '''
    this function is to caculate the similarity between x list and
    NRPS or PKS (Y).it is going to be used to compare detected cluster
    with NRPS and PKS.Since rearrange of genes within a gene cluster did
    happen in nature and did not change the function of the gene cluster.
    we believe the gene ordervof a gene cluster doesn't matter here.
    similarity between [0,1,2,3] and [0,1,2,4]is 0.75. similarity
    between [0,1,2,3,4] and [0,1,2,4] is 0.8.similarity between
    [0,1,2,4] and [0,1,2,3,4] is 1
    '''
    intersection_cardinality = len(set.intersection(*[set(x), set(y)]))
    return intersection_cardinality/float(len(set(x)))


def _remove_net_NRPS_PKS(hits, query_pfamnumbers):
    '''
    Because NRPS and PKS can be well detected by antismash, and net NRPS
    or PKS (not includeing addentional enzymes) are too frequently found
    in bacterial genomes and are not interesting for us, we remove them
    from the result. But if a cluster with NRPS or PKS together with
    addentional enzymes was shared by different genomes, it will be
    reported
    '''
    net_NRPS = {501, 13193, 550, 668, 13745, 501}
    net_PKS = {109, 550, 106, 107, 8240, 14765, 698, 2801, 108, 8990, 2797, 
               975, 13193, 501}
    # combine two sets if hybrid RRPS+PKS are also too often found 
    # and not interesting to us
    refined_hits = []
    for hit in hits:
        hit_pfamnumbers = query_pfamnumbers[hit.query_start: hit.query_end]
        if (_belong_similarity(hit_pfamnumbers, net_NRPS) < 0.8 
            and _belong_similarity(hit_pfamnumbers, net_PKS) < 0.8):
            refined_hits.append(hit)
    return refined_hits


def _filter_query_dna_length(query_g_pfam_list, hits, dna_length_threshold):
    """ 
    Removes all hits, where the dna sequence length of the query sublist
    is not long enough.
    """
    filtered_hits = list()
    for hit in hits:
        query_dna_start = query_g_pfam_list[hit.query_start].start
        query_dna_end = query_g_pfam_list[hit.query_end - 1].end
        if (query_dna_end - query_dna_start) > dna_length_threshold:
           filtered_hits.append(hit)
    return filtered_hits


def _find_hits(query_seq_rec, ref_seq_rec, hit_param):
    """
    """
    query_seq = query_seq_rec._get_pfam_num_seq()
    ref_seq = ref_seq_rec._get_pfam_num_seq()

    hits = _find_common_sublists(query_seq, ref_seq, hit_param.seed_size)
    hits = _completely_group_neighboring_hits(hits, hit_param.gap_threshold)
    hits = _filter_by_size(hits, hit_param.size_threshold)
    hits = _remove_net_NRPS_PKS(hits, query_seq)
    hits = _filter_query_dna_length(query_seq_rec.g_pfam_list, 
                                   hits, hit_param.dna_length_threshold)
    return hits


def _write_hits_to_db(conn, hits, query_seq_rec, ref_seq_rec):
    '''
    Inserts the two sublists of the hit into the database.
    '''
    query_g_pfam_list = query_seq_rec.g_pfam_list
    ref_g_pfam_list = ref_seq_rec.g_pfam_list

    hit_rows = list()
    for hit in hits:
        hit_rows.append([
            query_seq_rec.seq_acc, 
            query_g_pfam_list[hit.query_start].pfam_id,
            query_g_pfam_list[hit.query_end - 1].pfam_id,
            ref_seq_rec.seq_acc,
            ref_g_pfam_list[hit.ref_start].pfam_id,
            ref_g_pfam_list[hit.ref_end - 1].pfam_id
            ])

    with contextlib.closing(conn.cursor()) as c:
        sql = ''' INSERT INTO hits (query_seq_acc, 
                                    query_first_pfamID,
                                    query_last_pfamID, 
                                    ref_seq_acc, 
                                    ref_first_pfamID, 
                                    ref_last_pfamID) 
                       VALUES (?, ?, ?, ?, ?, ?) '''
        c.executemany(sql, hit_rows)


def _write_host_comparison_to_db(conn, query_hostID, ref_hostID):
    with conn:
        c = conn.cursor()
        sql = ''' INSERT INTO host_comparisons (query_hostID, ref_hostID) 
                       VALUES (?, ?) '''
        c.execute(sql,(query_hostID, ref_hostID))
        c.close()


def _find_hits_singlethread(host_pairs, hit_param, read_database,
                            write_database):
    """
    """
    # Connect to databases
    conn_read = sqlite3.connect(read_database)
    conn_write = sqlite3.connect(write_database)
    sqlite3.register_adapter(np.int64, lambda val: int(val))
    # HostID can be in numpy.int64 format. Sqlite3 can only read
    # numpy.int64 if it has been registered before
    
    previous_query_hostID = None
    previous_ref_hostID = None
    for query_hostID, ref_hostID in host_pairs:
        if _host_comparison_exists(conn_read, query_hostID, ref_hostID) \
            is False:

            # Get the seq records of each host
            if query_hostID != previous_query_hostID:
                query_seq_records = _get_seq_records(conn_read, query_hostID)
                previous_query_hostID = query_hostID
            if ref_hostID != previous_ref_hostID:
                ref_seq_records = _get_seq_records(conn_read, ref_hostID)
                previous_ref_hostID = ref_hostID
            # In many/most cases the query_hostID stays the same. We
            # do the same for the ref just in case the code changes in 
            # future

            # There might be more than one seq record per file
            seq_record_pairs = itertools.product(query_seq_records, 
                                                 ref_seq_records)
            # Compare the two hosts with each other
            with conn_write:
                for query_seq_rec, ref_seq_rec in seq_record_pairs:
                    hits = _find_hits(query_seq_rec, ref_seq_rec, hit_param)           
                    _write_hits_to_db(conn_write, hits, query_seq_rec, 
                                      ref_seq_rec)
                _write_host_comparison_to_db(conn_write, query_hostID, 
                                             ref_hostID)

    conn_read.close()
    conn_write.close()


def _combine_databases(main_database, temp_database):
    '''
    '''
    conn = sqlite3.connect(main_database)
    with conn:
        with contextlib.closing(conn.cursor()) as c:
            c = conn.cursor()
            c.execute(''' ATTACH DATABASE ? AS temp_db ''', (temp_database,))  
            c.execute('''
                INSERT OR IGNORE INTO hits (query_seq_acc, 
                                            query_first_pfamID, 
                                            query_last_pfamID, 
                                            ref_seq_acc, 
                                            ref_first_pfamID, 
                                            ref_last_pfamID)
                               SELECT query_seq_acc, 
                                      query_first_pfamID,
                                      query_last_pfamID,
                                      ref_seq_acc,
                                      ref_first_pfamID,
                                      ref_last_pfamID
                                 FROM temp_db.hits
                ''')  
            c.execute('''
                INSERT OR IGNORE INTO host_comparisons (query_hostID, 
                                                        ref_hostID)
                               SELECT query_hostID, 
                                      ref_hostID
                                 FROM temp_db.host_comparisons
                ''')
            c.execute(''' COMMIT ''')
            c.execute(''' DETACH DATABASE temp_db ''')  



def find_hits_multithread(host_pairs, hit_param, database_path, threads):
    """
    """
    try:
        _create_temporary_databases(database_path, threads)
        host_pair_chunks = _host_pairs_into_chunks(host_pairs, threads)
        Parallel(n_jobs=threads)(delayed(_find_hits_singlethread) \
            (host_pairs, hit_param, database_path, 
            _temp_database_path(database_path, thread)) \
                for thread, host_pairs in host_pair_chunks)
        for thread, _ in host_pair_chunks:
            temp_database_path = _temp_database_path(database_path, thread)
            _combine_databases(database_path, temp_database_path)
    
    except Exception:
        print_stderr(traceback.format_exc())
    
    finally:
        for thread in range(threads):
            temp_database_path = _temp_database_path(database_path, thread)
            os.remove(temp_database_path) 
