import contextlib

def create_database_tables(conn):
    '''
    Creates all necessary tables
    '''
    c = conn.cursor()
    c.execute('''
                CREATE TABLE "hosts" (
                	"hostID"	        INTEGER NOT NULL UNIQUE,
                	"organism"	        TEXT NOT NULL,
                	"host_type"	        TEXT,
                	"num_of_sequences"	INTEGER NOT NULL,
                	"L50"	            INTEGER NOT NULL,
                    "file"              TEXT NOT NULL UNIQUE,
                	PRIMARY KEY("hostID" AUTOINCREMENT)
                )
                '''
    )
    c.execute('''
                CREATE TABLE "seq_records" (
                    "seq_accession"	    TEXT NOT NULL UNIQUE,
					"hostID"	    TEXT NOT NULL,
	                "description"	TEXT NOT NULL,
	                "seq_length"	INTEGER NOT NULL,
	                FOREIGN KEY("hostID") REFERENCES "hosts"("hostID"),
	                PRIMARY KEY("seq_accession")
                )
                '''
    )

    c.execute('''
                CREATE TABLE "cds" (
                	"locus_tag"	    TEXT NOT NULL UNIQUE,
					"seq_accession"	TEXT NOT NULL,
                	"product"	    TEXT NOT NULL,
                	"translation"	TEXT NOT NULL,
                	"cds_start"	    INTEGER NOT NULL,
                	"cds_end"	    INTEGER NOT NULL,
                	"core_genome"   INTEGER,
                	FOREIGN KEY("seq_accession") REFERENCES \
                        "seq_records"("seq_accession"),
                	PRIMARY KEY("locus_tag")
                )
                '''
    )
    c.execute('''
                CREATE TABLE "pfams" (
                	"pfamID"	        INTEGER NOT NULL UNIQUE,
                	"pfam_num"	        TEXT NOT NULL,
                    "locus_tag"         TEXT NOT NULL,
					"seq_accession"	    TEXT NOT NULL,
                	"pfam_start"	    INTEGER NOT NULL,
                	"pfam_end"	        INTEGER NOT NULL,
                	"strand"	        INTEGER NOT NULL,
                    "antismash_core"    INTEGER NOT NULL,
                    "transporter"       INTEGER,
	                FOREIGN KEY("locus_tag") REFERENCES "cds"("locus_tag"),
                    FOREIGN KEY("seq_accession") REFERENCES \
                        "seq_records"("seq_accession"),
	                PRIMARY KEY("pfamID" AUTOINCREMENT)
                )
                '''
    )
    c.execute('''
                CREATE TABLE "sublists" (
                	"clusterID"	            INTEGER,
                	"host_type"	            TEXT NOT NULL,
                	"first_pfamID"	        INTEGER NOT NULL,
                	"last_pfamID"		    INTEGER NOT NULL,
                    "core_genome_fraction"  REAL NOT NULL,
                    "antismash_fraction"    REAL NOT NULL,

                	FOREIGN KEY("first_pfamID") REFERENCES "pfams"("pfamID"),
                    FOREIGN KEY("last_pfamID") REFERENCES "pfams"("pfamID"),
                    PRIMARY KEY("first_pfamID", "last_pfamID")
                )
                '''
    )
    c.execute('''
                CREATE TABLE "hits" (
                    "hitID"                 INTEGER NOT NULL UNIQUE,
                	"query_first_pfamID"	INTEGER NOT NULL,
                    "query_last_pfamID"		INTEGER NOT NULL,
                	"ref_first_pfamID"		INTEGER NOT NULL,
                    "ref_last_pfamID"		INTEGER NOT NULL,
                    PRIMARY KEY("hitID" AUTOINCREMENT)
                )
                '''
    )
    c.execute('''
                CREATE TABLE "cluster" (
                	"clusterID"	                INTEGER NOT NULL UNIQUE,
	                "max_core_genome_fraction"	REAL,
	                "transporter_indicator"	    REAL,
	                "number_core_pfams"	        INTEGER,
                    "min_antismash_fraction"    REAL,
                    "max_antismash_fraction"    REAL,
                	PRIMARY KEY("clusterID" AUTOINCREMENT)
                )
                '''
    )
    c.execute('''
                CREATE TABLE "host_comparisons" (
                	"query_hostID"	INTEGER NOT NULL,
                	"ref_hostID"	INTEGER NOT NULL,
                	FOREIGN KEY("query_hostID") REFERENCES "hosts"("hostID"),
                	FOREIGN KEY("ref_hostID") REFERENCES "hosts"("hostID"),
                	PRIMARY KEY("query_hostID", "ref_hostID")
                )
                '''
    )
    c.close()


def add_database_indeces(conn):
    """
    Improves the speed of data retrieval operations on a database table
    at the cost of additional writes and storage space to maintain the
    index data structure. Indexes are used to quickly locate data
    without having to search every row. 

    Maybe the indece creation should moved to the places where each 
    respective index is used.    
    """
    with contextlib.closing(conn.cursor()) as c:

        # Used in transp.py
        c.execute(''' 
                CREATE INDEX idx_pfams_pfamnum
                ON pfams (pfam_num) 
                ''')   

        # Used in compare.py
        c.execute(''' 
                CREATE INDEX idx_pfams_seqaccession
                ON pfams (seq_accession) 
                ''')
        c.execute(''' 
                CREATE INDEX idx_seqrecords_hostID
                ON seq_records (hostID)
                ''')
        c.execute('''
                CREATE INDEX idx_hosts_hosttype_L50
                ON hosts (host_type, L50)
                ''')

        # Used in cluster.py
        c.execute(''' 
                CREATE INDEX idx_hits_query_pfamIDs
                ON hits (query_first_pfamID, query_last_pfamID) 
                ''')
        # Used in cluster.py
        c.execute(''' 
                CREATE INDEX idx_hits_ref_pfamIDs
                ON hits (ref_first_pfamID, ref_last_pfamID) 
                ''')
        # Used in cluster.py and results.py
        c.execute(''' 
                CREATE INDEX idx_sublists_clusterID
                ON sublists (clusterID) 
                ''')
        
        # Used in results.py
        c.execute(''' 
                CREATE INDEX idx_cluster_indicators
                ON cluster (core_genome_indicator, 
                            antismash_indicator, 
                            transporter_indicator
                            ) 
                ''')
