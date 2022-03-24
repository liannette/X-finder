import os.path
import urllib.request
from goatools.base import get_godag
from goatools.gosubdag.gosubdag import GoSubDag
import re
from general_functions import get_base_dir

def get_transporter_pfams_from_pfam2go():
    '''
    Returns all pfams from pfam2go that are associated to a GO_ID which has the 
    molecular_function "transporter activity (GO:0005215) or the 
    biological_process "localization (GO:0051179)" as an "is a"/"part_of"
    ancestor
    '''
    # Initialize
    transport_pfams = set()
    # get GO graph
    godag = get_godag('go-basic.obo', optional_attrs='relationship')

    # get pfam2go file through url
    file = urllib.request.urlopen("http://current.geneontology.org/ontology/external2go/pfam2go")
    lines = file.readlines()
    for line in lines:
        line = line.decode("utf-8")
        # find lines that associate a pfam to a GO_ID
        if line.startswith("Pfam:"):
            line = line.split()
            pfam = int(re.search(r"PF(\d{5})", line[0]).group(1)) #int(line[0][7:])
            GO_ID = line[-1]
            #print(pfam, GO_ID)

            # Get all ancestors of the GO_ID
            optional_relationships = {'part_of'}    # "is_a" already included
            gosubdag_r1 = GoSubDag([GO_ID], godag, relationships=optional_relationships, prt=None)
            try:
                ancestors = gosubdag_r1.rcntobj.go2ancestors[GO_ID]
            except:
                #print("Unexpected error:", sys.exc_info()[0])
                continue
            #print(GO_ID)

            # Check if the GO ID is (part of) a localization process or transport function
            if "GO:0051179" in ancestors or "GO:0005215" in ancestors:
                transport_pfams.add(pfam)
                #print('{GO} ancestors: {P}'.format(
                #    GO=GO_ID,
                #    P=gosubdag_r1.rcntobj.go2ancestors[GO_ID]))
    transport_pfams = sorted(transport_pfams)
    return transport_pfams


def main():

    # Get transporter pfams
    transporter_pfams = get_transporter_pfams_from_pfam2go()

    filename = "transporter_pfams.txt"
    outfile_path = os.path.join(get_base_dir(), "data", "transporter pfams", "transporter_pfams.txt")
    with open(outfile_path, "w") as outfile:
        for pfam in transporter_pfams:
            print(pfam, file=outfile)


if __name__ == "__main__":
    main()