# Consider changing the code so that ref/query information is added 
# at a later stage.
# This would allow to have only one database containing all the genomes
# and the comparisons would then be saved in a seperate database.

import glob
from Bio import SeqIO
import os.path
import traceback
from xfinder.common import print_stdout, print_stderr
import sys
import contextlib
import re
import sqlite3


class Pfam:

    def __init__(self, pfam_num, cds_locus_tag, start, end, strand):
        self.pfam_num = pfam_num
        self.cds_locus_tag = cds_locus_tag
        self.start = start  # 0-based
        self.end = end      # 0-based
        self.strand = strand
        self.antismash_core = 0

    @classmethod
    def _from_gbk_feature(cls, feature):
        """ Used in importgbk.py """
        pfam_feature = cls(
            cls._get_pfam_number(feature.qualifiers["db_xref"][0]),
            feature.qualifiers["locus_tag"][0],
            feature.location.nofuzzy_start,
            feature.location.nofuzzy_end,
            feature.strand,
        )
        return pfam_feature

    @staticmethod
    def _get_pfam_number(pfam_id):
        pfam_num = re.search(r"PF(\d{5})(:?$|\D)", pfam_id)
        try:
            if pfam_num is None:
                raise Exception("Pfam number not found")
        except Exception:
            print_stderr("Failed to extract Pfam number from PfamID '{}'\n{}"
                         "".format(pfam_id, traceback.format_exc()))
        return int(pfam_num.group(1))

    def _pfam_overlap(self, previous_pfam):
        """ 
        Checks is the pfam overlaps with another pfam. The other pfam 
        needs to have a lower start value than self. Returns true, if 
        the overlap is bigger than 1/2 of either pfam size, otherwise
        returns False
        """
        if previous_pfam.start > self.start:
            raise ValueError(f'previous_pfam.start > pfam.start')
        overlap = previous_pfam.end - self.start
        if (overlap > (self.end - self.start)/2 or 
            overlap > (previous_pfam.end - previous_pfam.start)/2):
            return True
        else:
            return False

    def _to_list(self):
        return [self.cds_locus_tag, self.pfam_num, self.start, self.end,
                self.strand, self.antismash_core]


class Cds:

    def __init__(self, locus_tag, product, translation, start, end):
        self.locus_tag = locus_tag
        self.product = product
        self.translation = translation
        self.start = start  # 0-based
        self.end = end      # 0-based

    @classmethod
    def _from_gbk_feature(cls, feature):
        cds_feature = cls(
            feature.qualifiers["locus_tag"][0],
            feature.qualifiers["product"][0],
            feature.qualifiers["translation"][0],
            feature.location.nofuzzy_start,
            feature.location.nofuzzy_end,
        )
        return cds_feature

    def _to_list(self):
        return [self.locus_tag, self.product, self.translation, self.start, 
                self.end]


class ProtoCore:
    def __init__(self, start, end):
        self.start = start  # 0-based
        self.end = end      # 0-based

    @classmethod
    def _from_gbk_feature(cls, feature):
        proto_core = cls(
            feature.location.nofuzzy_start,
            feature.location.nofuzzy_end
        )
        return proto_core

    def _contains_pfam(self, pfam):
        if (self.start <= pfam.start 
            and self.end >= pfam.end):
            return True
        else:
            return False


class SeqRecord:

    def __init__(self, acc, organism, description, seq_length, cds_list, g_pfam_list):
        self.seq_acc = acc
        self.organism = organism
        self.description = description
        self.seq_length = seq_length
        self.cds_list = cds_list
        self.g_pfam_list = g_pfam_list # genome pfam sequence

    @classmethod
    def _from_gbk_feature(cls, rec):
        """
        accession (incl version), describtion and organism 
        are extracted from the sequence record; all CDS records and the
        genome pfam list are extracted. Not all pfams are in the 
        genome pfam list; if two pfams have a overlap area bigger then
        half of either of them, only the one with smallest pfam number 
        will be kept.
        """
        cds_list = list()
        g_pfam_list = list()
        antismash_core_list = list()

        for feature in rec.features:
            if feature.type == "proto_core":
                antismash_core_list.append(
                    ProtoCore._from_gbk_feature(feature)
                    )
            elif feature.type == "CDS":
                try:
                    cds = Cds._from_gbk_feature(feature)
                except KeyError:
                    # record does not contain all needed qualifiers
                    continue
                else:
                    cds_list.append(cds)
            elif feature.type == "PFAM_domain": 
                pfam = Pfam._from_gbk_feature(feature)
                # g_pfam_list is without overlapping (>50%) pfams
                cls._add_pfam_feature(g_pfam_list, pfam)

        for proto_core in antismash_core_list:
            for pfam in g_pfam_list:
                if proto_core._contains_pfam(pfam):
                    pfam.antismash_core = 1

        record = cls(
            acc=rec.id, # is ACCESSION.VERSION in gbk file
            organism=rec.annotations["organism"], 
            description=rec.description, # is DEFINITION in gbk file
            seq_length=len(rec.seq), 
            cds_list=cds_list, 
            g_pfam_list=g_pfam_list,
        )
        return record

    @staticmethod
    def _add_pfam_feature(g_pfam_list, pfam):
        """
        For pfam features overlapping with each other, only the one with
        smallest pfam number will be kept. (alternatively we can choose 
        the one with highest score (not done yet)). However, in case of
        overlaps smaller than 1/2, both pfam features will be tolerated.
        To tolerate gene inversion, the direction of pfam is not 
        concidered.
        """
        # Some PFAM_domain features don't have PFAM ID (this is a bug 
        # from antismash (they use a dictionary of pfam name and pfam ID
        # which need to be updated mannually), remember to check if this
        # happen too often, if yes ask Kai to fix it). 
        # ADD exceptions!
        if len(g_pfam_list) == 0:
            g_pfam_list.append(pfam)
        else:
            previous_pfam = g_pfam_list[-1]
            if pfam._pfam_overlap(previous_pfam):
                # If >=50% overlap with previous pfam, 
                # keep only the pfam with smallest pfam_num
                if pfam.pfam_num < previous_pfam.pfam_num: 
                    g_pfam_list.pop()
                    g_pfam_list.append(pfam)
            else:
                # no overlap, add pfam 
                g_pfam_list.append(pfam)


def list_of_gbk_files(genome_dirs):
    genome_files = list()
    for genome_dir in genome_dirs:
        genome_files += glob.glob(os.path.join(genome_dir, "*.gbk"))
    return genome_files


def _l50(records):
    '''
    Calculates the L50 value of a list of records. The L50 value
    indicates how many records (actually contigs) are needed to cover at
    least 1/2 of the combined sequence length.
    '''
    seq_lengths = sorted([rec.seq_length for rec in records], reverse=True)
    partial_seq_length = 0
    l50 = 0
    while partial_seq_length < sum(seq_lengths)/2:
        partial_seq_length += seq_lengths[l50]
        l50 += 1
    return l50


def _host_to_db(conn, records, host_type, gbk_file):
    """
    Inserts a new host into the table hosts, returns the hostID (rowid)
    """
    organism = records[0].organism
    file = os.path.basename(gbk_file)
    with contextlib.closing(conn.cursor()) as c:

        sql = ''' INSERT INTO hosts (organism, host_type, num_of_sequences, L50, 
                                     file) 
                       VALUES (?,?,?,?,?) '''
        c.execute(sql, (
            organism, host_type, len(records), _l50(records), file))
        hostID = c.execute('''SELECT last_insert_rowid()''').fetchall()[0][0]
        # hostID is rowid() column:
        # "If the table has a column of type INTEGER PRIMARY KEY then that
        # column is another alias for the rowid"
        # (http://www.sqlite.org/c3ref/last_insert_rowid.html)

    return hostID


def _records_to_db(conn, records, host_type, gbk_file):
    """
    Inserts the record into the database. It is assumed that all records
    are from the same host/organism, therefore only one entry is made in
    the host table.
    """
    with conn:
        # context manager ensures that any transaction are rolled back
        # if any exception occurs, or committed otherwise
        hostID = _host_to_db(conn, records, host_type, gbk_file)
        with contextlib.closing(conn.cursor()) as c:
            for rec in records:
                sql = ''' INSERT INTO seq_records (hostID, seq_accession, 
                                                   description, seq_length)
                               VALUES (?,?,?,?) '''
                c.execute(sql, (hostID, rec.seq_acc, rec.description, 
                                rec.seq_length))
                sql = ''' INSERT INTO cds (seq_accession, locus_tag, product,
                                        translation, cds_start, cds_end)
                               VALUES (?,?,?,?,?,?) '''
                c.executemany(sql, [[rec.seq_acc] + cds._to_list() 
                                    for cds in rec.cds_list])
                sql = ''' INSERT INTO pfams (seq_accession, locus_tag, 
                                             pfam_num, pfam_start, pfam_end, 
                                             strand, antismash_core) 
                               VALUES (?,?,?,?,?,?,?) '''
                c.executemany(sql, [[rec.seq_acc] + pfam._to_list() 
                                for pfam in rec.g_pfam_list])
        c.close()


def gbk_file_to_db(conn, gbk_file, host_type):
    """
    """
    records = [SeqRecord._from_gbk_feature(rec) 
               for rec in SeqIO.parse(gbk_file, "gb")]
    try:
        if len(set([rec.organism for rec in records])) > 1:
            # We expect all sequences to be from the same organism
            raise Exception("Different organisms in the same gbk file")
        _records_to_db(conn, records, host_type, gbk_file)
        return True
    except Exception:
        print_stderr("An error occured while importing the file '{}'. \n{}"
                     "".format(os.path.basename(gbk_file), 
                               traceback.format_exc()))
        return False


def print_import_progress(finished_files, total_files):
    """ 
    Returns True if the number of finished files is a multiple of 20 or 
    when the number of finished files equals the total number of files.
    Otherwise returns False.
    """
    if finished_files % 50 == 0 or finished_files == total_files:
        return True
    else:
        return False


def import_genomes(database_path, host_type, genome_dirs, out_dir):
    try:
        gbk_files = list_of_gbk_files(genome_dirs)
        print_stdout(f"Importing {len(gbk_files)} genbank files ({host_type}) "
                    "into the database.", out_dir)
        
        conn = sqlite3.connect(database_path)
        for i in range(len(gbk_files)):
            success = gbk_file_to_db(conn, gbk_files[i], host_type)
            if success is False:
                print_stdout("File was not imported into the database: "
                            f"{gbk_files[i]}", out_dir)
            if print_import_progress(i+1, len(gbk_files)):
                print_stdout("{}/{} files done".format(i+1, len(gbk_files)), 
                            out_dir)
        conn.close()

    except:
        print_stderr(traceback.format_exc(), out_dir)
        sys.exit(1)