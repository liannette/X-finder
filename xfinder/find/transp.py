import os
import tempfile
import re
import urllib.request
from goatools.base import get_godag
from goatools.gosubdag.gosubdag import GoSubDag
from xfinder.find.common import print_stdout
import contextlib
import sqlite3


def _get_transporter_pfams_from_pfam2go():
    '''
    Returns all pfams from pfam2go that are associated to a GO_ID which
    has the molecular_function "transporter activity (GO:0005215) or the 
    biological_process "localization (GO:0051179)" as an "is a" or
    "part_of" ancestor
    https://github.com/tanghaibao/goatools/blob/main/notebooks/parent_go_terms.ipynb
    '''
    transporter_pfams = set()

    with tempfile.TemporaryDirectory() as temp_dir:

        # Load the GO directed acyclic graph
        go_graph_path = os.path.join(temp_dir, 'go-basic.obo')
        go_graph = get_godag(go_graph_path, optional_attrs='relationship', 
                             prt=None) # will print some stuff anyway

        # This is the pfam2go file
        url = "http://current.geneontology.org/ontology/external2go/pfam2go"

        for line in urllib.request.urlopen(url).readlines():
            line = line.decode("utf-8")
            if line.startswith("Pfam:"):
                line = line.split()
                pfam_num = int(re.search(r"PF(\d{5})$", line[0]).group(1))
                go_id = re.search(r"(GO:\d{7})$", line[-1]).group(1)

                # Subgraph with "is_a" and "part_of" relationships
                # Returns an empty subgraph and prints a line if go_id
                # is not in go_graph. Main reason are obsolete terms.
                go_subgraph = GoSubDag([go_id], go_graph, 
                                       relationships={'part_of'}, prt=None)
                try:
                    ancestors = go_subgraph.rcntobj.go2ancestors[go_id]
                    if ("GO:0005215" in ancestors # transporter activity
                        or "GO:0051179" in ancestors): # localization
                        transporter_pfams.add(pfam_num)
                except KeyError:
                    # go_id was not found in empty go_subgraph
                    continue

    return transporter_pfams
    # print to file
    #with open(file_path, "w") as outfile:
    #    for pfam in sorted(transport_pfams):
    #        print(pfam, file=outfile)


def get_transporter_pfams(file_path=None):
    """ 
    Get all pfams associated with transporter activity 
    """
    if file_path is None:
        transporter_pfams = _get_transporter_pfams_from_pfam2go()
    elif os.path.exists(file_path) is False:
        print_stdout("File with transporter pfams does not exist yet. "
                     "Transporter pfams will be collected and written "
                     "to file.")
        transporter_pfams = _get_transporter_pfams_from_pfam2go()
    else:
        with open(file_path, "r") as infile:
            transporter_pfams = set(infile.read().splitlines())
    return transporter_pfams



def add_transporter_pfam_information(transporter_pfams, database_path):
    """ 
    Adds the transporter pfam information to the database table pfams. 
    """
    conn = sqlite3.connect(database_path)
    with conn:
        with contextlib.closing(conn.cursor()) as c:
            c.execute(''' UPDATE pfams SET transporter = 0 ''')
            sql = ''' UPDATE pfams 
                         SET transporter = 1 
                       WHERE pfam_num = ? '''
            c.executemany(sql, [(pfam,) for pfam in transporter_pfams])
    conn.close()
