import sqlite3
import os.path
import pandas as pd
from general_functions import connect_to_db, timepoint, showtimepoints, get_base_dir
import argparse
import sys
from datetime import datetime


# Functions --------------------------------------------------------------------


def index_db(conn):
    with conn:
        c = conn.cursor()
        c.execute(''' 
                CREATE INDEX cds_contig_idx
                ON cds (contig)
                ''')
        c.execute(''' 
                CREATE INDEX hits_clusterID_idx
                ON hits (clusterID)
                ''')
        c.close()


def drop_index(conn):
    with conn:
        c = conn.cursor()
        c.execute(''' DROP INDEX cds_contig_idx ''')
        c.execute(''' DROP INDEX hits_clusterID_idx ''')
        c.close()


def get_all_cluster(conn, core_genome_cutoff, transporter_cutoff):
    '''
    Get a list of all clusterIDs
    '''
    sql = '''   SELECT clusterID, core_genome_indicator, transporter_indicator, number_core_pfams
                FROM cluster 
                WHERE core_genome_indicator <= ?
                    AND transporter_indicator <= ?
                    AND  manual_exception IS NULL
                ''' 
    c = conn.cursor()
    c.execute(sql, (core_genome_cutoff, transporter_cutoff))
    return c.fetchall()


def get_sublists_in_cluster(conn, clusterID):
    '''
    '''
    c = conn.cursor()
    sql = """
            SELECT hits.query_pfamID_start, hits.query_pfamID_end, hits.ref_pfamID_start, hits.ref_pfamID_end
            FROM hits 
            WHERE clusterID = ?
          """
    rows = c.execute(sql, (clusterID,)).fetchall()
    q_start, q_end, r_start, r_end = zip(*rows)
    query_sublists = sorted(set(zip(["query"]*len(rows), q_start, q_end)))
    ref_sublists = sorted(set(zip(["ref"]*len(rows), r_start, r_end)))
    return query_sublists, ref_sublists


def get_locustags_pfamnums(conn, first_pfamID, last_pfamID):
    """
    """
    sql = """ 
            SELECT hosts.organism, hosts.file, contigs.contig, contigs.description, pfams.locus_tag, pfams.pfamnumber, pfams.pfamstart, pfams.pfamend
            FROM pfams
		        INNER JOIN contigs ON pfams.contig = contigs.contig
		        INNER JOIN hosts ON hosts.hostID = contigs.hostID
            WHERE pfams.pfamID BETWEEN ? AND ?
          """
    c = conn.cursor()
    rows = c.execute(sql, (first_pfamID, last_pfamID)).fetchall()
    organism, gbk_filename, seq_id, seq_description, cds_locustags, pfamnums, pfamstarts, pfamends = zip(*rows)
    sublist_start = min(pfamstarts)
    sublist_end = max(pfamends)
    #cds_locustags = sorted(set(cds_locustags)) # to get rid of duplicates, the original order is preserved because locus_tags are sequentially numbered
    first_last_cds_locustag = (cds_locustags[0], cds_locustags[-1])
    return organism[0], gbk_filename[0], seq_id[0], seq_description[0], sublist_start, sublist_end, first_last_cds_locustag, pfamnums


def get_cds_products(conn, first_last_cds_locustag, contig, underline):
        
    c = conn.cursor()

    # get the ROWID of the first and the last cds locustag
    sql = """
          SELECT ROWID
          FROM cds
          WHERE locus_tag = ?
          """ 
    first_cds_rowid = c.execute(sql, (first_last_cds_locustag[0],)).fetchall()[0][0]
    last_cds_rowid = c.execute(sql, (first_last_cds_locustag[1],)).fetchall()[0][0]
    
    #Get the cds of the sublist
    sql = """
          SELECT product, core_genome
          FROM cds
          WHERE ROWID BETWEEN {0} AND {1}
          """.format(first_cds_rowid, last_cds_rowid)
    rows = c.execute(sql).fetchall()
    cds_products, core_genome_labels = zip(*rows)

    # upstream and downstream cds
    sql = """
          SELECT product, core_genome
          FROM cds
          WHERE ROWID BETWEEN {0}-15 AND {0}-1
            AND contig = ?
          """.format(first_cds_rowid)
    rows = c.execute(sql, (contig,)).fetchall()
    if len(rows) == 0:
        us_cds_products = list()
        us_core_genome = list()
    else: 
        us_cds_products, us_core_genome = zip(*rows)
    us_cds_labelled = mark_core_genome_cds(us_cds_products, label_list=us_core_genome, method=underline)

    sql = """
          SELECT product, core_genome
          FROM cds
          WHERE ROWID BETWEEN {0}+1 AND {0}+15
            AND contig = ?
          """.format(last_cds_rowid)
    rows = c.execute(sql, (contig,)).fetchall()
    if len(rows) == 0:
        ds_cds_products = list()
        ds_core_genome = list()
    else:
        ds_cds_products, ds_core_genome = zip(*rows)
    ds_cds_labelled = mark_core_genome_cds(ds_cds_products, label_list=ds_core_genome, method=underline)

    return cds_products, core_genome_labels, us_cds_labelled, ds_cds_labelled


def create_cds_label_dict(cds_products, core_genome_labels):
    """
    Create a dict for all cds in cluster, key: cds_product, value: all occuring core genome labels 
    """
    cds_in_cluster_dict = dict()
    for i in range(len(cds_products)):           # Iterates through each sublist
        for j in range(len(cds_products[i])):    # Iterates through the cds
            cds_product = cds_products[i][j]
            core_genome_label = core_genome_labels[i][j]
            if cds_product in cds_in_cluster_dict.keys():
                cds_in_cluster_dict[cds_product].add(core_genome_label)
            else: 
                cds_in_cluster_dict[cds_product] = set(core_genome_label)
    return cds_in_cluster_dict


def mark_core_cluster_cds(cds_products, core_cds_in_cluster): 
    """
    Adds the core cluster cds label
    """  
    labelled_cds_products = list()
    for cds_product in cds_products:     
        if cds_product not in core_cds_in_cluster:
            labelled_cds_products.append(cds_product + "*")
        else: 
            labelled_cds_products.append(cds_product)
    return labelled_cds_products


def mark_core_genome_cds(cds_products, label_dict=None, label_list=None, method=None):
    """
    Adds the core genome label
    """
    labelled_cds_products = list()
    # For summary results
    if label_dict: 
        for cds_product in cds_products:
            label = "".join(sorted(label_dict[cds_product.rstrip("*")]))
            labelled_cds_products.append("{:>2} {}".format(label, cds_product)) 
    # for detailed results
    elif label_list:
        for i in range(len(cds_products)):
            if label_list[i] == "-":
                labelled_cds_products.append(method)
            labelled_cds_products.append(cds_products[i])
            labelled_cds_products.append(";\n") # add spacer
        labelled_cds_products.pop() # delete last spacer
    return labelled_cds_products


def print_summary_header(database, core_genome_cutoff, transporter_cutoff, outfile):
    """
    Prints a header for the summary outfile. Specifies the database name, the applied filters and labels.
    """
    header =    "Database: {}\n"\
                "Core genome indicator cutoff: {}\n" \
                "Transporter indicator cutoff: {}\n\n" \
                "Core genome indicator = highest fraction of core genome CDS in one sublist\n" \
                "Transporter indicator = fraction of core pfams (in cluster) associated with transporter function\n\n" \
                " + = does not align to the core genome of Streptomyces (<90% identity)\n" \
                " - = aligns to the core genome of Streptomyces (>90% identity)\n" \
                "+- = aligns in >=1 reference sublist and does not align in >= reference sublist \n" \
                "CDS product with neither + or - were only found in query sublists\n\n" \
                "* = CDS is not found in every sublists of the cluster".format(
                    database, core_genome_cutoff, transporter_cutoff
                )
    print_msg_box(msg=header, indent=2, outfile=outfile)


def print_msg_box(msg, outfile, indent=1, width=None, title=None):
    """
    Print message-box with optional title. 
    Taken from https://stackoverflow.com/questions/39969064/how-to-print-a-message-box-in-python/40080853
    """
    lines = msg.split('\n')
    space = " " * indent
    if not width:
        width = max(map(len, lines))
    box = f'╔{"═" * (width + indent * 2)}╗\n'  # upper_border
    if title:
        box += f'║{space}{title:<{width}}{space}║\n'  # title
        box += f'║{space}{"-" * len(title):<{width}}{space}║\n'  # underscore
    box += ''.join([f'║{space}{line:<{width}}{space}║\n' for line in lines])
    box += f'╚{"═" * (width + indent * 2)}╝'  # lower_border
    print(box, file=outfile)


def print_cluster_to_summary_results(entry, outfile): 
    '''
    Prints an cluster to the summary text file
    '''
    clusterID, core_genome_indicator, transporter_indicator, number_core_pfams, num_query_sublists, num_ref_sublists, cds_in_cluster = entry
    print("> Cluster {} \n" \
          "Number of sublists: {} ({} query, {} ref) \n" \
          "Core genome indicator: {} \n" \
          "Transporter indicator: {} \n" \
          "Number of core pfams: {}".format(clusterID, 
                                            num_query_sublists + num_ref_sublists, 
                                            num_query_sublists, 
                                            num_ref_sublists,
                                            core_genome_indicator, 
                                            transporter_indicator, 
                                            number_core_pfams), 
           file=outfile)
    print("CDS products:", file=outfile)
    for cds_product in cds_in_cluster:
        print('   ', cds_product, file=outfile)
    print("", file=outfile)


# Main program -----------------------------------------------------------------


def main(database, core_genome_cutoff, transporter_cutoff, outfile_prefix):
    
    detailed_results = list()
    summary_results = list()

    # Create  XlsxWriter Excel object
    detailed_results_path = os.path.join(get_base_dir(), "results", outfile_prefix + "_detailed.xlsx")
    writer = pd.ExcelWriter(detailed_results_path, engine='xlsxwriter')
    workbook  = writer.book
    underline = workbook.add_format({'underline': True})

    conn = connect_to_db(database)
    # Create index, for faster querying. skip if index already exists. This could be because the program was interruped previously
    try:
        index_db(conn)
    except:
        pass
    print("{}: Indexed database.".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))

    with conn:

        # Get a list of all cluster that satisfy the cutoff values
        cluster_list =  get_all_cluster(conn, core_genome_cutoff, transporter_cutoff)
        print("{}: Got a list of all cluster.".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
        

        for clusterID, core_genome_indicator, transporter_indicator, number_core_pfams in cluster_list:
            core_genome_indicator = round(core_genome_indicator, 3)
            transporter_indicator = round(transporter_indicator, 3)

            # Get all the sublists in the cluster
            query_sublists, ref_sublists = get_sublists_in_cluster(conn, clusterID)
            sublists = ref_sublists + query_sublists

            # Get all relevant information of the sublists
            sublists_detailed = list()
            for hosttype, first_pfamID, last_pfamID in sublists: 
                organism, gbk_filename, seq_id, seq_description, sublist_start, sublist_end, first_last_cds_locustag, pfamnums = get_locustags_pfamnums(conn, first_pfamID, last_pfamID)
                cds_products, core_genome_labels, us_cds_labelled, ds_cds_labelled = get_cds_products(conn, first_last_cds_locustag, seq_id, underline)
                sublists_detailed.append([
                    hosttype, 
                    gbk_filename,
                    organism, 
                    seq_id, 
                    seq_description, 
                    sublist_start, 
                    sublist_end, 
                    cds_products,
                    core_genome_labels,
                    pfamnums, 
                    us_cds_labelled, 
                    ds_cds_labelled])

            # Only the cds information
            cds_products = [sublist[7] for sublist in sublists_detailed]
            core_genome_labels = [sublist[8] for sublist in sublists_detailed]

            # Core cds in cluster
            core_cds_in_cluster = set(cds_products[0]).intersection(*cds_products)
            
            # For summary file: Create a dict for all cds in cluster, key: cds_product, value: all occuring core genome labels 
            cds_in_cluster_dict = create_cds_label_dict(cds_products, core_genome_labels)
            all_cds_in_cluster = sorted(cds_in_cluster_dict.keys(), key=lambda x: x.lower())
            all_cds_in_cluster = mark_core_cluster_cds(all_cds_in_cluster, core_cds_in_cluster)
            all_cds_in_cluster = mark_core_genome_cds(all_cds_in_cluster, label_dict=cds_in_cluster_dict)
            summary_results.append([clusterID, core_genome_indicator, transporter_indicator, number_core_pfams, len(query_sublists), len(ref_sublists), all_cds_in_cluster])

            for hosttype, gbk_filename, organism, seq_id, seq_description, sublist_start, sublist_end, cds_products, core_genome_labels, pfamnums, us_cds_labelled, ds_cds_labelled in sublists_detailed:
                cds_products_set = sorted(set(cds_products), key=lambda x: x.lower())
                cds_products_set = mark_core_cluster_cds(cds_products_set, core_cds_in_cluster)
                cds_products_set =  mark_core_genome_cds(cds_products_set, label_list=core_genome_labels, method=underline)
                cds_products = mark_core_cluster_cds(cds_products, core_cds_in_cluster)
                cds_products = mark_core_genome_cds(cds_products, label_list=core_genome_labels, method=underline)
                detailed_results.append(
                    [clusterID, len(sublists_detailed), core_genome_indicator, transporter_indicator, number_core_pfams,
                    hosttype, organism, gbk_filename, seq_id, seq_description, sublist_start, sublist_end, 
                    cds_products_set,  
                    cds_products, 
                    ", ".join(pfamnums), 
                    us_cds_labelled, 
                    ds_cds_labelled])

        print("{}: Got all sublists and the CDS products.".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))

    drop_index(conn)
    conn.close()
    print("{}: Dropped the database index again.".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))


    # Write detailed results to a tsv file
    detailed_results = pd.DataFrame(detailed_results, columns = [
        'clusterID', 'number of sublists', 'core genome indicator', 'transporter indicator', 'number of core pfams', 
        'hosttype', 'organism', 'file', 'sequence ID', 'sequence description', 'start position on contig', 'end position on contig', 
        'CDS products set', 'CDS products sequence', 'PFAM number sequence', 'upstream CDS products', 'downstream CDS products'
        ])
    detailed_results.to_excel(writer, sheet_name='Sheet1', index=False)
    
    worksheet = writer.sheets['Sheet1']
    worksheet.set_column(9, 9, 40)
    worksheet.set_column(12, 14, 60)
    worksheet.set_column(15, 17, 60)
    worksheet.autofilter(0, 0, len(detailed_results), len(detailed_results.columns) - 1)
    
    # underline core genome cds
    bg_color = "#F5F5F5"
    border_color = "#C0C0C0"
    cell_fill_format = workbook.add_format({'border':1, 'bg_color': bg_color, 'border_color':border_color})
    
    for row_idx in range(0, len(detailed_results)):

        if detailed_results.iloc[row_idx, 5] == "ref":
            worksheet.set_row(row_idx+1, 16, cell_fill_format)
            cds_cells_format = workbook.add_format({'border':1, 'text_wrap': True, 'border_color': border_color, 'bg_color': bg_color,})
        else:
            cds_cells_format = workbook.add_format({'border':1, 'text_wrap': True, 'border_color': border_color})

        for col_idx in [12, 13, 15, 16]:
            cds = detailed_results.iloc[row_idx, col_idx]
            if len(cds) >= 3: # write_rich_string() needs at least 3 fragments/formats
                worksheet.write_rich_string(row_idx+1, col_idx, *cds, cds_cells_format)
            elif len(cds) == 0: # no cds
                worksheet.write_blank(row_idx+1, col_idx, "", cds_cells_format)
            elif len(cds) == 1: # only one cds without underline
                worksheet.write_string(row_idx+1, col_idx, *cds, cds_cells_format)
            elif len(cds) == 2: # only one cds with underline
                cds_cells_format.set_underline()
                worksheet.write_string(row_idx+1, col_idx, cds[1], cds_cells_format)

    workbook.close()

    print("{}: Wrote the detailed results file".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))

    # Write summary results to txt file
    summary_results_path = os.path.join(get_base_dir(), "results", outfile_prefix + "_summary.txt")
    with open(summary_results_path, "w") as outfile_summary:
        print_summary_header(database, core_genome_cutoff, transporter_cutoff, outfile_summary)
        for entry in summary_results:
            print_cluster_to_summary_results(entry, outfile_summary)

    print("{}: Wrote the summary results file".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))




# Main program =================================================================

if __name__ == "__main__":

    print("Commandline input:", " ".join(sys.argv))

    parser = argparse.ArgumentParser()
    parser.add_argument('-db', action="store", dest="database", type=str, default="database.db", help='Name of the sql database (default=database.db')
    parser.add_argument('-core', action="store", dest="core_genome_cutoff", type=float, default=0.5, help='Cutoff fot the core genome indicator. From all ref sublists in the hits, is the max fraction of cds that aligned to the core genome. (default=0.5)')
    parser.add_argument('-transp', action="store", dest="transporter_cutoff", type=float, default=0.2, help='Cutoff for the transporter indicatior. Value indicates the fraction of transporter pfams in the core pfams of the cluster. (default=0.2)')
    parser.add_argument('-o', action="store", dest="outfile_prefix", type=str, default="results", help='Outfile suffix for the summary and detailed results. (default=results)')
    
    args = parser.parse_args()
    database = args.database
    core_genome_cutoff = args.core_genome_cutoff
    transporter_cutoff = args.transporter_cutoff
    outfile_prefix = args.outfile_prefix

    print("{}: Program started. Database: {}, core_genome_cutoff: {}, transporter_cutoff: {}, outfile_prefix: {}" \
        .format(
            datetime.now().strftime("%d/%m/%Y %H:%M:%S"), 
            database,
            core_genome_cutoff,
            transporter_cutoff,
            outfile_prefix
            )
        )
    main(database=database, core_genome_cutoff=core_genome_cutoff, transporter_cutoff=transporter_cutoff, outfile_prefix=outfile_prefix)
    print("{}: Program finished".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))