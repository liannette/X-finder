import os
import git
import datetime
import sys
import traceback
import re
import contextlib


class Pfam:

    def __init__(self, pfam_num, cds_locus_tag, start, end, strand, 
                 pfam_id=None):
        self.pfam_id = pfam_id
        self.pfam_num = pfam_num
        self.cds_locus_tag = cds_locus_tag
        self.start = start  # 0-based
        self.end = end      # 0-based
        self.strand = strand

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
                self.strand]


class Cds():

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
                

class SeqRecord:

    def __init__(self, acc, organism=None, description=None, seq_length=None, 
                 antismash_core_list=None, cds_list=None, g_pfam_list=None):
        self.seq_acc = acc
        self.organism = organism
        self.description = description
        self.seq_length = seq_length
        self.antismash_core_list = antismash_core_list
        self.cds_list = cds_list
        self.g_pfam_list = g_pfam_list # genome pfam sequence

    @classmethod
    def _from_db(cls, seq_acc, conn):
        ''' Returns the genome pfam sequence of a sequence record '''
        with contextlib.closing(conn.cursor()) as c:
            sql = ''' SELECT pfam_num, locus_tag, pfam_start, pfam_end, strand,
                             pfamID
                        FROM pfams 
                       WHERE seq_acc=? '''
            rows = c.execute(sql,(seq_acc,)).fetchall()
        return cls(acc=seq_acc, 
                   g_pfam_list=[Pfam(*row) for row in rows])

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
                start = feature.location.nofuzzy_start
                end = feature.location.nofuzzy_end
                antismash_core_list.append((start, end))
            elif feature.type == "CDS":
                cds_list.append(Cds._from_gbk_feature(feature))
            elif feature.type == "PFAM_domain":
                # g_pfam_list without overlapping (>50%) pfams 
                pfam = Pfam._from_gbk_feature(feature)
                cls._add_pfam_feature(g_pfam_list, pfam)

        record = cls(
            acc=rec.id, # is ACCESSION.VERSION in gbk file
            organism=rec.annotations["organism"], 
            description=rec.description, # is DEFINITION in gbk file
            seq_length=len(rec.seq), 
            antismash_core_list=antismash_core_list,
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

    def _get_pfam_num_seq(self):
        return [pfam.pfam_num for pfam in self.g_pfam_list]


def get_git_root():
    ''' Returns the absolute path of the repository root directory '''
    git_repo = git.Repo(__file__, search_parent_directories=True)
    git_root = git_repo.git.rev_parse("--show-toplevel")
    return git_root


def get_database_dir():
    ''' Returns the absolute path of the database directory '''
    database_dir = os.path.join(get_git_root(), "data", "database")
    return database_dir


def _now():
    return datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")


def print_stdout(msg):
    ''' Adds a timestamp to a string and prints it to stdout'''
    print("{}\t{}".format(_now(), msg))


def print_stderr(msg):
    ''' Adds a timestamp to a string and prints it to stderr'''
    print("{}\t{}".format(_now(), msg), file=sys.stderr)
