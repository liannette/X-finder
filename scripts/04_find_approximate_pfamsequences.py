#clearned up with the help of Simon Shaw <sisha@biosustain.dtu.dk>
# compared to version 1, 1)here we added a function to write the results into a genbank file. the results are the predicted gene clusters, they are added as features into the genbank file which usually is the file of the query genome. reference https://caretdashcaret.com/2016/06/15/how-to-create-genbank-files-with-biopython/ https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/genbank/
#2) we also updated the core function. so "one fragment from either of the two sequence can have more than one matches on the other, and all of them will be reported." instead of "if there are more than two matches in reference from one fragment of query, only the longest/first one will be kept.if two fragments of the query match one area in reference, both will be kept."
# here we assume all genomes are in linear form. this may overlook the genes near the breaking point of circular genomes.it should not be a problem if circular genomes are always break at the ori area. because the genes in ori area are already well studied.
# there are hard coding here. remember to change them when necessary
#the order of filters is important, try to change it when necessary.
#pfam database has updated to deal with overlapped pfam annotation. try to find out more.

import sqlite3
import sys
import os.path
import numpy as np
from general_functions import connect_to_db, timepoint, showtimepoints
from joblib import Parallel, delayed
from datetime import datetime
import argparse
import itertools
import pandas as pd


# Classes ======================================================================

class Hit:
    def __init__(self, query_start, query_end, ref_start, ref_end, hit_direction):
        assert query_start < query_end
        # since it is inside a class, when will it be run?
        self.query_start = query_start
        self.query_end = query_end
        assert ref_start < ref_end
        self.ref_start = ref_start
        self.ref_end = ref_end
        self.hit_direction = hit_direction  # direct == 1, inverted == -1

    def get_query_length(self):
        return self.query_end - self.query_start

    def get_ref_length(self):
        return self.ref_end - self.ref_start

    def has_same_query_location(self, other):
        return self.query_start == other.query_start and self.query_end == other.query_end

    def is_inside_query(self, other):
        return self.query_start >= other.query_start and self.query_end <= other.query_end

    def query_overlaps(self, other):
        return self.query_start < other.query_end and self.query_end > other.query_start

    def ref_overlaps(self, other):
        return self.ref_start < other.ref_end and self.ref_end > other.ref_start

    def both_regions_overlap(self, other):
        return self.query_overlaps(other) and self.ref_overlaps(other)

    def gap_calculation(self, other):
        '''calculate the gap between two sublists (sublist is called as hit in function find_common_sublists)
        a sublist/hit is represented by a tuple in the form of
        (start in query, end in query/not included, start in reference, end in reference/not included, direction*length).
        a gap value below 0 indicate that the two sublist overlap with eachother
        '''
            # here we calculate the distance/gap between the two sublists.
            #gaps within the two sublists were not considered
        gap_in_query = max(other.query_start-self.query_end, self.query_start-other.query_end)
        gap_in_reference = max(other.ref_start-self.ref_end, self.ref_start-other.ref_end)
        return gap_in_query, gap_in_reference
        # retrun a tuple
    def combine_sublist(self, other):
        '''
        combine two sublist into one.
        a sublist/hit is represented by a tuple in the form of
        (start in query, end in query/not included, start in reference, end in reference/not included, direction*length).
        the direction will be the sum of the directions of all the matched pfams
        don't try to add gap value in to the sublist tuple. when two sublist overlap with eachother, there is no way to calculate
        the gap value in the combined_sublist, as that require the position of the gaps.
        '''
        combined_start_in_query = min(other.query_start, self.query_start)
        combined_end_in_query = max(other.query_end, self.query_end)
        combined_start_in_reference = min(other.ref_start, self.ref_start)
        combined_end_in_reference = max(other.ref_end, self.ref_end)
        combined_direction = other.hit_direction + self.hit_direction
        combined = Hit(combined_start_in_query, combined_end_in_query, \
        combined_start_in_reference, combined_end_in_reference, combined_direction)
        return combined

    def __repr__(self):
        return str((self.query_start, self.query_end, self.ref_start, self.ref_end, self.hit_direction))


# Functions ====================================================================

def host_comparison_exists(conn, query_host, ref_host, times):
    '''
    Checks if hits have already been calculated for the two hits
    '''
    start_time = datetime.now()
    c = conn.cursor()
    sql = ''' SELECT COUNT(*) 
                FROM host_comparisons 
               WHERE query_host=? 
                 AND ref_host=?
          '''
    result = c.execute(sql,(query_host, ref_host)).fetchall()[0][0]
    c.close()
    times[2] += (datetime.now() - start_time).total_seconds()
    return result == 1, times

def read_contigs_from_db(conn, host, times):
    '''
    read the contigs belonging to a specific host from table contigs
    :param conn:
    :param host:
    :return: contigs
    '''
    start_time = datetime.now()
    sql = ''' SELECT contig 
                FROM contigs 
                WHERE host=?'''
    c = conn.cursor()
    rows = c.execute(sql,(host,)).fetchall()
    contigs = [row[0] for row in rows]
    c.close()
    times[3] += (datetime.now() - start_time).total_seconds()
    return contigs, times

def read_gpfamsequence_from_db(conn, contig, times):
    '''
    read a gpfamsequence from table pfams as comparison refernce
    :param conn:
    :param contig:
    :return: a python list of table rows
    '''
    start_time = datetime.now()
    sql = ''' SELECT ID, pfamnumber, locus_tag, pfamstart, pfamend 
                FROM pfams 
               WHERE contig=? '''
    c = conn.cursor()
    gpfamsequence = c.execute(sql,(contig,)).fetchall()
    c.close()
    times[4] += (datetime.now() - start_time).total_seconds()
    return gpfamsequence, times

def find_common_sublists(s1, s2, length_threshold, times):
    '''
    purpose: compare a query gpfamnumbersequence with a reference gpfamnumbersequence to find out all
     common pfamnumber sublists longer than the threshold (1 by default) in both directions.
    the common sublists are all exact matches. insertion, deletion or invertion was not tolerated in this step.
    these represent the exactly matched pfam domain sequences from both genomes.
    pfamnumber repeats (for example from pks and nrps genes) generate overlapped common_sublists.they will be combined in the next step.
    input: a list of pfamnumbers representing the query genome and a list represent the reference genomes. a length_threshold of the match, it should be 1 or 2.
    output: a list of hits. each hit is represented by a tuple in the form of (start in query, end in query/not included,
    start in reference, end in reference/not included, direction*length)
    a single match has a direction of 0, a match with length of 2 has a direction of 1 or -1
    potential bug. The direction of mach of substrings with tandom repeat of one element can be wrong.
    a match with length of 1
    '''
    start_time = datetime.now()
    #print("length of query genome pfam list:", len(s1))
    #print("length of reference genome pfam list:", len(s2))
    #print (s1[0])
    sublist_length_threshold = length_threshold
    #single pfam number match is also considered. test sublist_length_threshold = 2
    #s1 = list("jianglinpdclin")
    #s2 = list("gnaijlincpoggljl")
    #s1 = [101,"a","b",201,"a","b","a",102,103,104,"a","b",202,"a","b","a",105,106,107,108,109,110]
    #s2 = [1,2,3,"b","a","b",4,5,6,"b","a","b",7,8,9,10]
    #extract pfamnumber_list from gpfamsequence
    # substring_length_threshold by default it is 2. don't change it. we will do filtering later by length of substring or others.
    score = np.zeros((2+len(s1), 2+len(s2)), dtype=int)
    #put query s1 on the Y axis of the 2D arrary; put reference S2 on the X axis of the 2D arrary
    score_i = np.zeros((2+len(s1), 2+len(s2)), dtype=int)
    # "+2", because we need an initional value 0 on both side of the strings
    # score_i is for detection of inverted match, there must be two arrarys! to avoid
    hits = []
    # hit was represented by a tuple (siq,eiq,sir,eir,di), standing for start in query, end in query(not include), start in reference, end in reference, direct hit(1) or inverted hit (-1).
    for x in range(1, 1 + len(s1)):
        '''here we use the classic Dynamic programming for finding common substrings'''
        for y in range(1, 1 + len(s2)):
            if s1[x - 1] == s2[y - 1]:
                score[x][y] = score[x - 1][y - 1] + 1
                score_i[x][y] = score_i[x - 1][y + 1] + 1

                if (x == len(s1) or y == len(s2) or s1[x] != s2[y]) and score[x][y] >= sublist_length_threshold:

                    # if ( it reaches the right or bottom of the array or continued match stops) and the match length >= the threshold
                    # don't check s1[x] != s2[y] before x == len(s1) or y == len(s2) or else you will have a list index error
                    # because when x reach len(S1), S1[X] does not exist
                    #print "substring_length: %d" % score[x][y]
                    #print "substring in s1:", s1[x-score[x][y]:x]
                    #print "substring in s2:", s2[y-score[x][y]:y]
                    hits.append(Hit(x-score[x][y], x, y-score[x][y], y, score[x][y]-1))
                    #print "-" * 10
                if (x == len(s1) or y == 1 or s1[x] != s2[y-2]) and score_i[x][y] >= sublist_length_threshold :
                    # if (continued match stops or it reaches the left or bottom of the array) and the match length >= the threshold
                    #print "substring_length: %d" % score_i[x][y]
                    #print "substring in s1:", s1[x-score_i[x][y]:x]
                    #print "substring in s2:", s2[y-1:y+score_i[x][y]-1]
                    hits.append(Hit(x-score_i[x][y], x, y-1, y+score_i[x][y]-1, -1*(score_i[x][y]-1)))
                    #print "-" * 10
                else:
                    continue
                    # the hits list has a order of a increasing Y
    times[5] += (datetime.now() - start_time).total_seconds()
    return hits, times


def group_neighboring_sublists(sublists, gap_threshold):
    # https://en.wikipedia.org/wiki/Approximate_string_matching
    '''purpose: 1) to combine the sublists that are close to (including overlap with) eachother
    on BOTH the qurery and reference genomes. this is similar to but more complex than the classic "merge overlapping intervals found in a list" problem
    2)to group the sublists that are close to (including overlap with) eachother
    on both the qurery and reference genomes.
    one such group will represent a gene cluster with some degree of modification, such as insertion
    deletion or invertion.
    a group is a list of sublists.
    Overlapped sublists generated from pfamnumber repeats
    (for example from pks and nrps genes which have tandom repeats of pfam domains) will be kept in the group.
    The one pfam to one pfam corresponding relationship in a group maybe used in futher for generating figures of gene cluster alignment.
    However a group can be too complex to be interpreted by human eyes.
    3)gap_threshold is the biggest gap allowed.'''
    combined_sublists = []
    '''groups = []'''
    # groups is a list of groups
    # a group is a list of neighboring sublists
    # a sublists is a tuple representing (start in query, end in query/not included, start in reference, end in reference/not included, direction)
    to_skip = [0]
    #this is to avoid repeated combining.
    # for example if we have combined sublists[5], sublists[7] and sublists[8]. we don't want to combine 5 and 7 again.
    # we are doing the combining backwards
    #any element that can be combined with [7] must have been combined with retangle [7+8]
    #any element that can be combined with [5] must have been combined with retangle [5+7+8], so we should not start from [7] or [5] again.
    #[6] may be combined wtih something else. so it will still be checked.
    # we find a situation that [8] can not combine with [6] but after [8] combined with [7] and [5], [6] can combine with [5++78]. we solved this by function
    # completely_group_neighboring_sublists(sublists, gap_threshold)
    for  i in range(1,len(sublists)+1):
        # here we start from the last element, we compare it with all the element before it with a reversed order.
        if i in to_skip:
            # if the last element in this round has been combined before, we skip it.
            continue
        else:
            combined_sublist = sublists[-i]
            #use the last element in this round as the default/initial value
            '''group = [sublists[-i],]'''
            for  n in range(1,len(sublists)+1-i):
                gap_in_query, gap_in_reference = combined_sublist.gap_calculation(sublists[-i-n])
                # you can do this from a tuple. gap_calculation(combined_sublist, sublists[-i-n]) returns a tuple
                if gap_in_query > gap_threshold:
                    break
                        # sublists is orderred with increasing sublist[1] values (end in query/not included).
                        # so once sublists[-i-n] is too far away on qurery, all the following elements will be the same
                elif gap_in_reference <= gap_threshold:
                    # the two elements are not far away on BOTH query and reference
                    combined_sublist = combined_sublist.combine_sublist(sublists[-i-n])
                    #print("the sublist obsorbed:", sublists[-i-n])
                    #print("combined sublist:    ", combined_sublist)
                    '''group.append(sublists[-i-n])'''
                    to_skip.append(i+n)
                    # add this new obsobed element to the to_skip list
                else:
                    continue

            combined_sublists.append(combined_sublist)
            '''group.reverse()
            groups.append(group)'''

    combined_sublists.reverse()
    #print ("to skip list:")
    #print (to_skip)
    return combined_sublists


def completely_group_neighboring_sublists(sublists, gap_threshold, times):
    """
    we find a situation that [8] can not combine with [6] but after [8] combined with [7] and [5],
    [6] can combine with [5++78]. we solved this by doing the grouping repeatly till there is no more to group
    """
    start_time = datetime.now()
    # a list of the sublist tuple
    #grouping_round = 0
    for i in range(1,len(sublists)):
        #grouping_round = grouping_round+1
        listlength = len(sublists)
        sublists = group_neighboring_sublists(sublists, gap_threshold)
        #print (list)
        if len(sublists) == listlength:
            # no grouping happened in this round
            break
    #print ("After", grouping_round, "rounds of grouping" )
    times[6] += (datetime.now() - start_time).total_seconds()
    return sublists, times


def filter_by_size(hits, size_threshold, times):
    start_time = datetime.now()
    result_list = []
    for hit in hits:
        if hit.get_query_length() >= size_threshold \
        and hit.get_ref_length() >= size_threshold:
        #combined_sublist is not too short on both of the genomes
        # and the gap accumulation is not too long on both genomes
            result_list.append(hit)
    times[7] += (datetime.now() - start_time).total_seconds()
    return result_list, times


def remove_net_NRPS_PKS(hits, query_pfamnumbers, times):
    '''Because NRPS and PKS can be well detected by antismash, and net NRPS or PKS
    (not includeing addentional enzymes) are too frequently found in bacterial genomes
    and are not interesting for us, we remove them from the result. But if a cluster
    with NRPS or PKS together with addentional enzymes was shared by different genomes, it will be reported'''
    def belong_similarity(x,y):
        '''this function is to caculate the similarity between x list and NRPS or PKS (Y).
        it is going to be used to compare detected cluster with NRPS and PKS.
         Since rearrange of genes within a gene cluster did happen in nature and
         did not change the function of the gene cluster. we believe the gene order
         of a gene cluster doesn't matter here. similarity between [0,1,2,3] and [0,1,2,4]
         is 0.75. similarity between [0,1,2,3,4] and [0,1,2,4] is 0.8.
         similarity between [0,1,2,4] and [0,1,2,3,4] is 1 '''
        intersection_cardinality = len(set.intersection(*[set(x), set(y)]))
        return intersection_cardinality/float(len(set(x)))
    start_time = datetime.now()
    net_NRPS = {501, 13193, 550, 668, 13745, 501}
    net_PKS = {109, 550, 106, 107, 8240, 14765, 698, 2801, 108, 8990, 2797, 975, 13193, 501}
    # if hybrid RRPS+PKS are also too often found and not interesting to us. combine this two sets.
    refined0_hits = []
    for hit in hits:
        hit_pfamnumbers = query_pfamnumbers[hit.query_start: hit.query_end]
        #print ("hit.query_start: hit.query_end worked")
        if belong_similarity(hit_pfamnumbers, net_NRPS) < 0.8 and belong_similarity(hit_pfamnumbers, net_PKS) < 0.8:
            refined0_hits.append(hit)
    times[8] += (datetime.now() - start_time).total_seconds()
    return refined0_hits, times


def remove_hit_of_too_short_DNAlength(query, hits, DNA_length_threshold, times):
    start_time = datetime.now()
    result = [hit for hit in hits if (query[(hit.query_end-1)][4]-query[hit.query_start][3]) > DNA_length_threshold]
    times[9] += (datetime.now() - start_time).total_seconds()
    return result, times
        #hit.query_end-1 is the last pfam of the hit!!!!
        # rewrite pfams into class like what we did with hit.


def insert_hits_to_db(conn, query_contig, ref_contig, query_gpfamsequence, reference_gpfamsequence, hits, times):
    '''
    Inserts the two sublists of the hit into the database. Checks if a sublist
    already exists in the database and assignes them to clusters accordingly.
    '''
    start_time = datetime.now()
    with conn:
        c = conn.cursor()
        c.execute('BEGIN TRANSACTION') 
        # Put multiple database writes into one transacion. Instead of writing to (and locking) the file each and 
        # every time a write query is issued, the write will only happen once when the transaction completes.
        # "Unless already in a transaction, each SQL statement has a new transaction started for it. This is very 
        # expensive, since it requires reopening, writing to, and closing the journal file for each statement. 
        # This can be avoided by wrapping sequences of SQL statements with BEGIN TRANSACTION; and END TRANSACTION; statements. 
        # This speedup is also obtained for statements which don't alter the database. 
        # https://web.archive.org/web/20140331055857/http://web.utk.edu/~jplyon/sqlite/SQLite_optimization_FAQ.html

        # Write the hits to the table "hits"
        for hit in hits:
            #pfamIDs_q = [str(item[0]) for item in query_gpfamsequence[hit.query_start:hit.query_end]]
            #pfamIDs_r = [str(item[0]) for item in reference_gpfamsequence[hit.ref_start:hit.ref_end]]
            sql = ''' INSERT INTO hits (query_contig, query_pfamID_start, query_pfamID_end, ref_contig, ref_pfamID_start, ref_pfamID_end) 
                                     VALUES (?, ?, ?, ?, ?, ?) '''
            c.execute(sql, (
                query_contig, 
                str(query_gpfamsequence[hit.query_start][0]), 
                str(query_gpfamsequence[(hit.query_end-1)][0]),
                ref_contig, 
                str(reference_gpfamsequence[hit.ref_start][0]),
                str(reference_gpfamsequence[(hit.ref_end-1)][0])
                ))
        # End transaction and close cursor. The keyword COMMIT is a synonym for END TRANSACTION.
        c.execute('COMMIT')
        c.close()
    times[10] += (datetime.now() - start_time).total_seconds()
    return times


def insert_host_comparisons_to_DB(conn, query_host, ref_host):
    '''
    Writes the two hosts that have been compared to the database
    '''
    with conn:
        c = conn.cursor()
        # Write the two hosts of the comparison to the table "host_comparisons"
        sql = ''' INSERT INTO host_comparisons (query_host, ref_host) 
                                 VALUES (?, ?) '''
        c.execute(sql,(query_host, ref_host))
        c.close()


def find_hits(query_host, ref_hosts, database, seed_size, gap_threshold, size_threshold, DNA_length_threshold, thread_no):
    """
    """

    conn = connect_to_db(database)
    conn2 = connect_to_db(database.split('.')[0] + '_temp' + str(thread_no) + '.db')

    # Register np.int64 type, so that sqlite can read them as int
    sqlite3.register_adapter(np.int64, lambda val: int(val))

    times = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

    query_contigs, times = read_contigs_from_db(conn, query_host, times)
    query_gpfamsequence_contigs = list()
    query_pfamnumbers_contigs = list()
    for query_contig in query_contigs: 
        query_gpfamsequence, times = read_gpfamsequence_from_db(conn, query_contig, times)
        query_gpfamsequence_contigs.append(query_gpfamsequence)
        query_pfamnumbers_contigs.append([item[1] for item in query_gpfamsequence])

    for ref_host in ref_hosts:

        comparison_exists, times = host_comparison_exists(conn, query_host, ref_host, times)
        if comparison_exists is False:
            #print("{}: Comparison between query '{}' and ref '{}' does not exist in the database yet.".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), query_host, ref_host))

            ref_contigs, times = read_contigs_from_db(conn, ref_host, times)
            
            for query_idx, ref_idx in itertools.product(range(len(query_contigs)), range(len(ref_contigs))):


                query_contig = query_contigs[query_idx]
                query_gpfamsequence = query_gpfamsequence_contigs[query_idx] # = list of pfam entries (pfamID, pfamnumber, locus_tag, pfamstart, pfamend)
                query_pfamnumbers = query_pfamnumbers_contigs[query_idx]
                ref_contig = ref_contigs[ref_idx]
                reference_gpfamsequence, times = read_gpfamsequence_from_db(conn, ref_contig, times) # = list of pfam entries (pfamID, pfamnumber, locus_tag, pfamstart, pfamend)
                ref_pfamnumbers = [item[1] for item in reference_gpfamsequence]
                #print("{}: Got the pfams for both hosts.".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S)))

                # Find the hits
                exact_matches, times = find_common_sublists(query_pfamnumbers, ref_pfamnumbers, seed_size, times)
                #print("{}: find_common_sublists finished.".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
                approximate_matches, times = completely_group_neighboring_sublists(exact_matches, gap_threshold, times)
                #print("{}: completely_group_neighboring_sublists finished.".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
                large_approximate_matches, times = filter_by_size(approximate_matches, size_threshold, times)
                #print("{}: filter_by_size finished.".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
                hits, times = remove_net_NRPS_PKS(large_approximate_matches, query_pfamnumbers, times)
                #print("{}: remove_net_NRPS_PKS finished.".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
                hits, times = remove_hit_of_too_short_DNAlength(query_gpfamsequence, hits, DNA_length_threshold, times)
                #print("{}: remove_hit_of_too_short_DNAlength finished. Will write results to db now.".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
                
                # Write hits/matches to database
                times = insert_hits_to_db(conn2, query_contig, ref_contig, query_gpfamsequence, reference_gpfamsequence, hits, times)
            
            insert_host_comparisons_to_DB(conn2, query_host, ref_host)
        
    conn.close()
    conn2.close()
    return times


def count_comparisons(database):
    '''
    Counts how many entries there are in the table host_comparisons
    '''
    conn = connect_to_db(database)
    sql = ''' SELECT COUNT(*) FROM host_comparisons '''
    c = conn.cursor()
    c.execute(sql,)
    count = c.fetchall()[0][0]
    c.close()
    conn.close()
    return count


def get_hosts(database, hosttype, max_L50):
    '''
    '''
    conn = connect_to_db(database)
    with conn:
        c = conn.cursor()
        sql = ''' SELECT host
                  FROM hosts 
                  WHERE hosttype = ? AND L50 <= ?'''
        rows = c.execute(sql, (hosttype, max_L50)).fetchall()
        c.close() 
    conn.close()
    return [row[0] for row in rows]


def create_temporary_databases(database, threads):
    '''
    '''
    for thread in range(threads):
        temp_database = database.split('.')[0] + '_temp' + str(thread) + '.db'
        conn = connect_to_db(temp_database)
        c = conn.cursor()
        c.execute('''
            CREATE TABLE "hits" (
            	"query_contig"	        TEXT NOT NULL,
            	"query_pfamID_start"    INTEGER NOT NULL,
            	"query_pfamID_end"	    INTEGER NOT NULL,
            	"ref_contig"	        TEXT NOT NULL,
            	"ref_pfamID_start"	INTEGER NOT NULL,
            	"ref_pfamID_end"	INTEGER NOT NULL,
                PRIMARY KEY("query_pfamID_start","ref_pfamID_start")
                )
                '''
        )
        c.execute('''
                CREATE TABLE "host_comparisons" (
                	"query_host"	TEXT NOT NULL,
                	"ref_host"	    TEXT NOT NULL,
                	PRIMARY KEY("query_host","ref_host")
                )
                '''
        )            
        conn.close()


def combine_databases(database, threads, BASE_DIR):
    '''
    '''
    conn = connect_to_db(database)
    for i in range(threads):
        temp_database = os.path.join(BASE_DIR, "data", "database", database.split('.')[0] + '_temp' + str(i) + '.db')
        with conn:
            c = conn.cursor()
            c.execute('ATTACH DATABASE ? AS temp_db', (temp_database,))  
            c.execute('''
                        INSERT OR IGNORE INTO hits (query_contig, query_pfamID_start, query_pfamID_end, ref_contig, ref_pfamID_start, ref_pfamID_end)
                        SELECT query_contig, query_pfamID_start, query_pfamID_end, ref_contig, ref_pfamID_start, ref_pfamID_end
                        FROM temp_db.hits
                        '''
                        )  
            c.execute('''
                        INSERT OR IGNORE INTO host_comparisons (query_host, ref_host)
                        SELECT query_host, ref_host
                        FROM temp_db.host_comparisons
                        '''
                        )  
            c.execute('COMMIT')
            c.execute('DETACH DATABASE temp_db')  
            c.close()
        os.remove(temp_database) 


# Run code functions ------------------------------------------------------------


def main(seed_size, gap_threshold, size_threshold, DNA_length_threshold, threads, database):

    # Get base directory
    dirname = os.path.dirname
    BASE_DIR = dirname(dirname(os.path.abspath(__file__)))
    
    # Get the number of comparisons already in the database
    comparison_count_start = count_comparisons(database)

    # Get query and reference hosts
    query_hosts = get_hosts(database, "query", max_L50=3)
    ref_hosts = get_hosts(database, "reference", max_L50=3)
    print("{}: Collected {} query hosts and {} reference hosts, which amounts to {} host comparisions.".format(
        datetime.now().strftime("%d/%m/%Y %H:%M:%S"), 
        len(query_hosts), 
        len(ref_hosts), 
        len(query_hosts)*len(ref_hosts)
        ))

    # Create temporary databases to enable concurrent database writes. They will be merged later.
    create_temporary_databases(database, threads)
    print("{}: Created temporary databases to allow concurrent writes. Starting to find and write hits now.".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))

    # Split the reference hosts into n=threads chunks and enumerate
    ref_hosts_split = list(enumerate(np.array_split(ref_hosts, threads)))

    matching_times = list()
    now = datetime.now()
    start = now
    for i in range(len(query_hosts)):
        # Using multithreading to find hits
        times = Parallel(n_jobs=threads)(delayed(find_hits)(query_hosts[i], ref_hosts, database, seed_size, gap_threshold, size_threshold, DNA_length_threshold, thread_no) for thread_no, ref_hosts in ref_hosts_split)
        total_duration = datetime.now() - start
        print("{}: {}/{} batches done. Batch duration: {}, Total Duration: {}".format(
            datetime.now().strftime("%d/%m/%Y %H:%M:%S"), i+1, len(query_hosts), datetime.now() - now, total_duration))
        matching_times += times
        now = datetime.now()
    
    # Combine the temporary tables
    combine_databases(database, threads, BASE_DIR)
    # Calculate how many new comparions were made 
    comparison_count_new = count_comparisons(database) - comparison_count_start
    print("{}: Combined the results from the temporary databases.".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
    

    print("{}: Programm finished. {} new comparisons have been added to the database.".format(
        datetime.now().strftime("%d/%m/%Y %H:%M:%S"), 
        comparison_count_new
        ))


    # Create the pandas DataFrame
    df = pd.DataFrame(matching_times, columns = [
        "query", 
        "reference",
        "host_comparison_exists",
        "read_contigs_from_db",
        "get_gpfamsequences",
        "find_common_sublists",
        "completely_group_neighboring_sublists", 
        "filter_by_size", 
        "remove_net_NRPS_PKS",
        "remove_hit_of_too_short_DNAlength",
        "insert_hits_to_db"
        ])
    #df.to_csv('matching_times.csv', index=False) 

    print("\nTimes, in percent (*threads):")
    print(df.iloc[:, 2:].sum(axis = 0)/total_duration.total_seconds()*100)
            



# Main program =================================================================

if __name__ == "__main__":

    print("Commandline input:", " ".join(sys.argv))

    parser = argparse.ArgumentParser()
    parser.add_argument('-t', action="store", dest="threads", type=int, default=1, help='Number of threads for multiprocessing (default=1)')
    parser.add_argument('-db', action="store", dest="database", type=str, default="database.db", help='Name of the sql database (default=database.db')
    args = parser.parse_args()
    threads = args.threads
    database = args.database

    print("{}: Program started. Threads: {}, database: {}".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), threads, database))

    main(seed_size=2, 
         gap_threshold=2, 
         size_threshold=6, 
         DNA_length_threshold=7000, 
         threads=threads, 
         database=database
         )
    # to check if the mode is the main program. if it is the if condition is true then the main() will be excuted
    # if it is mode imported by other mode or the main program. __name__= the name of this file.
    #the if condition will be faulse, and the main() will not be excuted
    # https://stackoverflow.com/questions/419163/what-does-if-name-main-do


    # Regarding multithreading: SQLite accepts several connections at the same time, but only one connection can write at the same time.