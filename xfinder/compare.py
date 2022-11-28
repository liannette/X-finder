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
import numpy as np
from xfinder.common import print_stdout, print_stderr
import sys
from joblib import Parallel, delayed
import itertools
import traceback


class Pfam:

    def __init__(self, pfam_id, pfam_num, cds_locus_tag, start, end, strand):
        self.pfam_id = pfam_id
        self.pfam_num = pfam_num
        self.cds_locus_tag = cds_locus_tag
        self.start = start  # 0-based
        self.end = end      # 0-based
        self.strand = strand


class SeqRecord:

    def __init__(self, acc, g_pfam_list):
        self.seq_acc = acc
        self.g_pfam_list = g_pfam_list # genome pfam sequence

    @classmethod
    def _from_db(cls, seq_acc, conn):
        # Get the genome pfam sequence of a sequence record 
        with contextlib.closing(conn.cursor()) as c:
            sql = ''' 
                SELECT pfamID, pfam_num, locus_tag, pfam_start, pfam_end,
                       strand
                  FROM pfams 
                 WHERE seq_accession=? 
                 '''
            rows = c.execute(sql,(seq_acc,)).fetchall()
        return cls(acc=seq_acc, 
                   g_pfam_list=[Pfam(*row) for row in rows])


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

    def _to_db_rows(self, query_seq_rec, ref_seq_rec):
        query_row = [
            "query", 
            query_seq_rec.g_pfam_list[self.query_start].pfam_id,
            query_seq_rec.g_pfam_list[self.query_end - 1].pfam_id,
            self.query_antismash_ind,
            self.query_core_genome_ind,
        ]
        ref_row = [
            "ref",
            ref_seq_rec.g_pfam_list[self.ref_start].pfam_id,
            ref_seq_rec.g_pfam_list[self.ref_end - 1].pfam_id,
            self.ref_antismash_ind,
            self.ref_core_genome_ind,
        ]
        return (query_row, ref_row)

    def __repr__(self):
        rep = ("Hit(" 
               + ", ".join([self.query_start, self.query_end, self.ref_start,
                            self.ref_end, self.hit_direction]) 
               + ")")
        return rep


def count_comparisons(database_path):
    '''
    Counts how many entries there are in the table host_comparisons
    '''
    conn = sqlite3.connect(database_path)
    with contextlib.closing(conn.cursor()) as c: # auto-closes cursor
        sql = ''' SELECT COUNT(*) FROM host_comparisons '''
        c.execute(sql,)
        num_comparisons = c.execute(sql).fetchall()[0][0]
    conn.close()
    return num_comparisons


def get_host_ids(conn, host_type, max_l50):
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
        sql = ''' SELECT seq_accession 
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
    query_pfamseq = [pfam.pfam_num for pfam in query_seq_rec.g_pfam_list] 
    ref_pfamseq = [pfam.pfam_num for pfam in ref_seq_rec.g_pfam_list]

    hits = _find_common_sublists(query_pfamseq, ref_pfamseq, 
                                 hit_param.seed_size)
    hits = _completely_group_neighboring_hits(hits, hit_param.gap_threshold)
    hits = _filter_by_size(hits, hit_param.size_threshold)
    hits = _remove_net_NRPS_PKS(hits, query_pfamseq)
    hits = _filter_query_dna_length(query_seq_rec.g_pfam_list, 
                                   hits, hit_param.dna_length_threshold)
    return hits


def _add_filter_information(hits, query_gpfam_list, ref_gpfam_list, conn):
    # Fraction of pfams that were marked antismash core
    sql_antismash = '''
        SELECT AVG(antismash_core)
          FROM pfams
         WHERE pfamID >= ? AND pfamID <= ?
        '''
    # Fraction of cds that were marked as core genome
    sql_coregenome = '''
        WITH distict_cds AS (
            SELECT DISTINCT cds.locus_tag, core_genome
              FROM pfams
            	   INNER JOIN cds ON pfams.locus_tag=cds.locus_tag
             WHERE pfamID >= ? AND pfamID <= ?
            )
        SELECT AVG(core_genome)
          FROM distict_cds
         ''' 
    with contextlib.closing(conn.cursor()) as c:
        for hit in hits:
            # For the query
            query_first_pfamID = query_gpfam_list[hit.query_start].pfam_id
            query_last_pfamID = query_gpfam_list[hit.query_end - 1].pfam_id
            
            c.execute(sql_antismash, (query_first_pfamID, query_last_pfamID))
            antismash_ind = c.fetchall()[0][0]
            hit.query_antismash_ind = round(antismash_ind, 3)
            
            c.execute(sql_coregenome, (query_first_pfamID, query_last_pfamID))
            core_genome_ind = c.fetchall()[0][0]
            hit.query_core_genome_ind = round(core_genome_ind, 3)

            # for the reference
            ref_first_pfamID = ref_gpfam_list[hit.ref_start].pfam_id
            ref_last_pfamID = ref_gpfam_list[hit.ref_end - 1].pfam_id
            
            c.execute(sql_antismash, (ref_first_pfamID, ref_last_pfamID))
            antismash_ind = c.fetchall()[0][0]
            hit.ref_antismash_ind = round(antismash_ind, 3)
            
            c.execute(sql_coregenome, (ref_first_pfamID, ref_last_pfamID))
            core_genome_ind = c.fetchall()[0][0]
            hit.ref_core_genome_ind = round(core_genome_ind, 3)
            
    return hits


def _find_hits_singlethread(query_hostID, query_seq_records, ref_hostID, 
                            hit_param, database_path):
    # Connect to databases
    conn = sqlite3.connect(database_path)
    # HostID can be in numpy.int64 format. Sqlite3 can only read
    # numpy.int64 if it has been registered before
    sqlite3.register_adapter(np.int64, lambda val: int(val))
    
    if _host_comparison_exists(conn, query_hostID, ref_hostID) is False:

        # Get the seq records of the ref host
        ref_seq_records = _get_seq_records(conn, ref_hostID)
        
        # Compare the two hosts with each other
        all_hits = []
        for q_seq_rec, r_seq_rec in itertools.product(query_seq_records, 
                                                      ref_seq_records):
            hits = _find_hits(q_seq_rec, r_seq_rec, hit_param)
            hits = _add_filter_information(hits,
                                           q_seq_rec.g_pfam_list, 
                                           r_seq_rec.g_pfam_list,
                                           conn)
            all_hits.append((q_seq_rec, r_seq_rec, hits))

        conn.close()
        return {"hosts": (query_hostID, ref_hostID), "all_hits": all_hits}
    else:
        conn.close()


def _write_hits_to_db(conn, hits, query_seq_rec, ref_seq_rec):
    '''
    Inserts the two sublists of the hit into the database.
    '''
    with contextlib.closing(conn.cursor()) as c:
        for hit in hits:
            query_row, ref_row = hit._to_db_rows(query_seq_rec, ref_seq_rec)
            sql = '''
                INSERT OR IGNORE INTO sublists (host_type, 
                                                first_pfamID, 
                                                last_pfamID, 
                                                antismash_fraction,
                                                core_genome_fraction)
                               VALUES (?, ?, ?, ?, ?) 
                '''
            c.execute(sql, query_row)
            c.execute(sql, ref_row)

            sql = ''' 
                INSERT INTO hits (query_first_pfamID,
                                  query_last_pfamID, 
                                  ref_first_pfamID, 
                                  ref_last_pfamID) 
                    VALUES (?, ?, ?, ?) 
                '''
            c.execute(sql, query_row[1:3] + ref_row[1:3])


def _write_host_comparison_to_db(conn, query_hostID, ref_hostID):
    c = conn.cursor()
    sql = ''' INSERT INTO host_comparisons (query_hostID, ref_hostID) 
                   VALUES (?, ?) '''
    c.execute(sql,(query_hostID, ref_hostID))
    c.close()
    

def find_hits_multithread(query_hostIDs, ref_hostIDs, hit_param, database_path,
                          threads):
    conn = sqlite3.connect(database_path)
    # HostID can be in numpy.int64 format. Sqlite3 can only read
    # numpy.int64 if it has been registered before
    sqlite3.register_adapter(np.int64, lambda val: int(val))

    for query_hostID in query_hostIDs:
        query_seq_records = _get_seq_records(conn, query_hostID)
        
        # split ref hosts into chunks, so mem usage wont go crazy
        chunksize = 20
        for i in range(0, len(ref_hostIDs), chunksize):
            ref_hostIDs_chunk = ref_hostIDs[i:i + chunksize]

            results = Parallel(n_jobs=threads)(
                delayed(_find_hits_singlethread)(
                    query_hostID,
                    query_seq_records, 
                    ref_hostID, 
                    hit_param, 
                    database_path
                    ) for ref_hostID in ref_hostIDs_chunk
                )

            # For each thread
            for result in results:
                # If the comparison already existed, result is None
                if result is not None:
                    # roll back automatically if an sql error occurs
                    with conn:
                        query_hostID, ref_hostID = result["hosts"]
                        _write_host_comparison_to_db(
                            conn, 
                            query_hostID, 
                            ref_hostID
                            )
                        # For each combination of sequence records
                        for q_sr, r_sr, hits in result["all_hits"]:
                            _write_hits_to_db(
                                conn, 
                                hits, 
                                q_sr, 
                                r_sr
                                )
        
        
def compare_genomes(seed_size, gap_threshold, size_threshold, 
                    DNA_length_threshold, threads, max_l50, database_path, 
                    out_dir):
    try:
        # Get the number of comparisons already in the database
        num_comparisons_start = count_comparisons(database_path)
        
        # Get the hostIDs
        conn = sqlite3.connect(database_path)
        query_hostIDs = get_host_ids(
            conn, 
            "query", 
            max_l50
            )
        ref_hostIDs = get_host_ids(
            conn, 
            "ref", 
            max_l50
            )
        conn.close()
        
        print_stdout(
            f"Comparing {len(query_hostIDs)} query hosts and {len(ref_hostIDs)} "
            "reference hosts, which amounts to  "
            f"{len(query_hostIDs)*len(ref_hostIDs)} comparisions.", out_dir)


        # Put the hit parameter in a nice class
        hit_parameters = HitParameters(
            seed_size, 
            gap_threshold, 
            size_threshold, 
            DNA_length_threshold
            )

        # find hits using multithreading
        find_hits_multithread(
            query_hostIDs, 
            ref_hostIDs,
            hit_parameters, 
            database_path, 
            threads
            )
            
        # Calculate how many new comparions were made 
        comparison_count_new = count_comparisons(database_path) \
                            - num_comparisons_start
        print_stdout(f"{comparison_count_new} new comparisons have been added to "
                    "the database.", out_dir)
    except:
        print_stderr(traceback.format_exc(), out_dir)
        sys.exit(1)