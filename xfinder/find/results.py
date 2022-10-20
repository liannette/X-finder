import sqlite3
import os.path
import pandas as pd
import contextlib

# functions still need to be put into order

class Sublist:
    
    def __init__(self, sublistID, host_type, file, organism, seq_acc,
                 description, start, end):
        self.id = sublistID
        self.host_type = host_type
        self.file = file
        self.organism = organism
        self.seq_acc = seq_acc
        self.description = description
        self.start = start
        self.end = end

    @classmethod
    def _from_row(cls, row):
        sublist = cls(*row)
        return sublist
        
        
class Cds:
    
    def __init__(self, locus_tag, product, core_genome):
        self.locus_tag = locus_tag
        self.product = product
        self.core_genome = core_genome
        self.pfams = None
    
    @classmethod
    def _from_row(cls, row):
        cds = cls(*row)
        return cds
    

def get_filtered_cluster(core_genome_cutoff, antismash_cutoff, 
                         transporter_cutoff, database_path):
    conn = sqlite3.connect(database_path)
    with contextlib.closing(conn.cursor()) as c:
        sql = '''   
            SELECT clusterID, 
                   core_genome_indicator, 
                   antismash_indicator, 
                   transporter_indicator, 
                   number_core_pfams
              FROM cluster 
             WHERE core_genome_indicator <= ?
                   AND antismash_indicator <= ?
                   AND transporter_indicator <= ?
          ORDER BY clusterID
            ''' 
        c.execute(sql, (core_genome_cutoff, 
                        antismash_cutoff, 
                        transporter_cutoff))
        rows = c.fetchall()
    conn.close()
    return rows


def _get_sublists(conn, clusterID):
    with contextlib.closing(conn.cursor()) as c:
        sql = """
            SELECT sublists.ROWID AS sublistID,
                   h.host_type, 
                   h.file,
                   h.organism,
                   s.seq_accession, 
                   s.description,
                   MIN(p.pfam_start),
                   MAX(p.pfam_end)
              FROM sublists
                   INNER JOIN pfams AS p
                       ON pfamID BETWEEN first_pfamID AND last_pfamID
                   INNER JOIN seq_records AS s
                       ON p.seq_accession = s.seq_accession
                   INNER JOIN hosts AS h 
                       ON s.hostID = h.hostID
             WHERE clusterID = ?
          GROUP BY sublistID
		  ORDER BY h.host_type
            """
        rows = c.execute(sql, (clusterID,)).fetchall()
    sublists = [Sublist._from_row(row) for row in rows]
    # Get the cds and pfam information for each sublist
    for sublist in sublists: 
        sublist_cds, upstream_cds, downstream_cds = _get_cds(sublist, conn)
        sublist.all_cds = sublist_cds
        sublist.upstream = upstream_cds
        sublist.downstream = downstream_cds
    # core cds products are present in all sublists of a cluster
    _add_core_cds_products(sublists) 
    return sublists


def _get_cds_from_rowids(seq_acc, first_cds_rowid, last_cds_rowid, conn):
    with contextlib.closing(conn.cursor()) as c:
        sql = """
            SELECT locus_tag, product, core_genome
            FROM cds
            WHERE ROWID BETWEEN ? AND ?
              AND seq_accession = ?
            """.format(first_cds_rowid, last_cds_rowid)
        c.execute(sql, (first_cds_rowid, last_cds_rowid, seq_acc))
        rows = c.fetchall()
    return [Cds._from_row(row) for row in rows]


def _get_cds(sublist, conn):
    
    # Get a pandas df that connects each pfam num in sublist with its cds
    sql = '''
        SELECT c.ROWID as cds_rowid,
        	   c.locus_tag, 
	    	   p.pfam_num, 
	    	   p.antismash_core,
	    	   p.pfam_start,
	    	   p.pfam_end
          FROM sublists
               INNER JOIN pfams AS p
	    			ON pfamID BETWEEN first_pfamID AND last_pfamID
	    	   INNER JOIN cds AS c
	    			ON p.locus_tag = c.locus_tag
         WHERE sublists.ROWID = ?
         '''
    df = pd.read_sql_query(sql, conn, params=(sublist.id,))
    # Combine the pfam num with the antismash core binary
    df["pfams"] = list(zip(df["pfam_num"], df["antismash_core"]))
    
    first_cds_rowid = min(df["cds_rowid"])
    last_cds_rowid = max(df["cds_rowid"])

    # Get all cds in sublist (cds with no pfams are not in df)
    sublist_cds = _get_cds_from_rowids(
        sublist.seq_acc, 
        first_cds_rowid, 
        last_cds_rowid, 
        conn)
    # Now, add the pfams to each cds
    for cds in sublist_cds:
        # list of (pfam_num, antismash_core_binary)
        cds.pfams = list(df[df["locus_tag"]==cds.locus_tag]["pfams"])

    # Get upstream and downstream cds
    upstream_cds = _get_cds_from_rowids(
        sublist.seq_acc, 
        first_cds_rowid - 15, 
        first_cds_rowid - 1, 
        conn)
    downstream_cds = _get_cds_from_rowids(
        sublist.seq_acc, 
        last_cds_rowid + 1, 
        last_cds_rowid + 15, 
        conn)
    
    return sublist_cds, upstream_cds, downstream_cds
    

def _add_core_cds_products(sublists):
    # core cds products are present in all sublists
    cds_products = list() # [cds.product for sublist in sublists for cds in sublist.all_cds]
    for sublist in sublists:
        cds_products.append([cds.product for cds in sublist.all_cds])
    core_cds_products = set(cds_products[0]).intersection(*cds_products)
    
    for sublist in sublists:
        for cds in sublist.all_cds:
            if cds.product in core_cds_products:
                cds.core_cds_product = 1
            else:
                cds.core_cds_product = 0
        for cds in sublist.upstream:
            if cds.product in core_cds_products:
                cds.core_cds_product = 1
            else:
                cds.core_cds_product = 0
        for cds in sublist.downstream:
            if cds.product in core_cds_products:
                cds.core_cds_product = 1
            else:
                cds.core_cds_product = 0


def _num_sublist_of_host_type(sublists, host_type):
    cnt = 0
    for sublist in sublists:
        if sublist.host_type == host_type:
            cnt += 1
    return cnt


def _core_genome_dict(sublists):
    """
    Create a dict for all cds in cluster, key: cds_product, value: all 
    occuring core genome labels (0 or 1)
    """
    core_genome_dict = dict()
    for sublist in sublists:
        for cds in sublist.all_cds:
            if cds.product in core_genome_dict.keys():
                core_genome_dict[cds.product].add(cds.core_genome)
            else: 
                core_genome_dict[cds.product] = set([cds.core_genome])
    return core_genome_dict


def _mark_cds_summary(cds_list, core_genome_dict):
    labeled1 = list()
    for cds_product, core_cds_product in cds_list:
        if core_cds_product == 0:
            # cds product not found in every sublist of cluster
            labeled1.append(cds_product + "*")
        else:
            labeled1.append(cds_product)
            
    labeled2 = list()
    for cds_product in labeled1:
        core_genome = core_genome_dict[cds_product.rstrip("*")]
        # some cds with this product align to core genome, some not
        if 0 in core_genome and 1 in core_genome:
            label = "+-"
        # all cds with this product align to core genome
        elif 1 in core_genome:
            label = "-"
        # no cds with this product align to core genome
        elif 0 in core_genome:
            label = "+"
        labeled2.append("{:>2} {}".format(label, cds_product)) 
    return labeled2
        

def _labeled_summary_cds_products(sublists):
    
    summary_cds = set()
    for sublist in sublists:
        for cds in sublist.all_cds:
            summary_cds.add((cds.product, cds.core_cds_product))
    sorted_summary_cds = sorted(summary_cds, key=lambda x: x[0].lower())
    core_genome_dict = _core_genome_dict(sublists)
    labeled_summary_cds = _mark_cds_summary(sorted_summary_cds, 
                                           core_genome_dict)
    return labeled_summary_cds
        

def _labeled_cds_products_set(sublist):
    
    cds_set = set()
    for cds in sublist.all_cds:
        cds_set.add((cds.product, cds.core_cds_product))
    sorted_cds_set = sorted(cds_set, key=lambda x: x[0].lower())
    
    # Only core cds products are labeled. 
    # To label core genome is not possible, because we have a set here,
    # and there might be two cds with the same product, one core genome
    # and one not. For the detailed results we have only the underline 
    # method to mark core genome cds.
    labeled_cds_set = list()
    for cds_product, core_cds_product in sorted_cds_set:
        if core_cds_product == 0:
            # cds product not found in every sublist of cluster
            labeled_cds_set.append(cds_product + "*")
        else:
            labeled_cds_set.append(cds_product)
         # add spacer
        labeled_cds_set.append(";\n")
    # delete last spacer
    labeled_cds_set.pop()
    
    return labeled_cds_set      


def _labeled_cds_products(sublist, core_genome_format, antismash_format):
    labelled_cds_products = list()
    for cds in sublist.all_cds:
        # Add core genome label
        if cds.core_genome == 1:
            # cds aligned to core genome
            labelled_cds_products.append(core_genome_format)
        # Add core cds product label
        if cds.core_cds_product == 0:
            # cds product not found in every sublist of cluster
            labelled_cds_products.append(cds.product + "*")
        else:
            labelled_cds_products.append(cds.product) 
            
        # Add pfams
        for i in range(len(cds.pfams)):
            pfam_num, antismash = cds.pfams[i]
            if i == 0:
                labelled_cds_products.append(" (")
            if antismash == 1:
                # pfam is in antismash core cluster
                labelled_cds_products.append(antismash_format)
            labelled_cds_products.append(pfam_num)
            labelled_cds_products.append(", ")
        # Delete last spacer of the pfam num
        if len(cds.pfams) > 0:
            labelled_cds_products.pop()
            labelled_cds_products.append(")")
        
        # add spacer for next cds
        labelled_cds_products.append(";\n")
        
    # delete last spacer of cds
    labelled_cds_products.pop()
    
    return labelled_cds_products


def _labeled_updownstream(sublist, core_genome_format, stream):
    
    if stream == "upstream":
        cds_list = sublist.upstream
    elif stream == "downstream":
        cds_list = sublist.downstream
        
    labelled_cds_products = list()
    for cds in cds_list:
        # Add core genome label
        if cds.core_genome == 1:
            # cds aligned to core genome
            labelled_cds_products.append(core_genome_format)
        labelled_cds_products.append(cds.product) 
        # add spacer for next cds
        labelled_cds_products.append(";\n")
    # delete last spacer of cds
    if len(labelled_cds_products) > 0:
        labelled_cds_products.pop()
    
    return labelled_cds_products


def _print_msg_box(msg, outfile, indent=1, width=None, title=None):
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


def _print_summary_header(database, core_genome_cutoff, antismash_cutoff, 
                          transporter_cutoff, outfile):
    """
    Prints a header for the summary outfile. Specifies the database name, the applied filters and labels.
    """
    #"CDS product with neither + or - were only found in query sublists\n\n" \
    header =    "Database: {}\n"\
                "Core genome indicator cutoff: {}\n" \
                "Antismash indicator cutoff: {}\n" \
                "Transporter indicator cutoff: {}\n\n" \
                "Core genome indicator = for each sublist, the fraction of CDS that align to the core genome of \n" \
                "                        Streptomyces is calculated. The core genome indicator is the highest \n" \
                "                        core genome fraction of all sublists in the cluster.\n" \
                "Antismash indicator   = for each sublist, the fraction of pfams that are in the core region of  \n" \
                "                        an antismash result is calculated. The antismash indicator is the highest \n" \
                "                        antismash fraction of all sublists in the cluster.\n" \
                "Transporter indicator = the fraction of core pfams (core pfams = pfam numbers that are present in \n" \
                "                        every sublist of the cluster) that are associated with transporter function\n\n" \
                " + = the CDS for this product does not align to the core genome of Streptomyces (<90% identity)\n" \
                " - = the CDS for this product aligns to the core genome of Streptomyces (>90% identity)\n" \
                "+- = at least one CDS for this product aligns to the core genome of Streptomyces, \n" \
                "     but also at least one CDS for this product aligns to the core genome of Streptomyces \n" \
                " * = CDS product is not found in every sublists of the cluster".format(
                    database, 
                    core_genome_cutoff, 
                    antismash_cutoff, 
                    transporter_cutoff,

                )
    _print_msg_box(msg=header, indent=2, outfile=outfile)


def _print_cluster_to_summary_results(entry, outfile): 
    '''
    Prints an cluster to the summary text file
    '''
      
    clusterID = entry[0]
    core_genome_indicator = entry[1]
    antismash_indicator = entry[2]
    transporter_indicator = entry[3]
    number_core_pfams = entry[4]
    num_query_sublists = entry[5]
    num_ref_sublists = entry[6]
    cds_in_cluster = entry[7]
            
    print("> Cluster {} \n" \
          "Number of sublists: {} ({} query, {} ref) \n" \
          "Core genome indicator: {} \n" \
          "Antismash indicator: {} \n" \
          "Transporter indicator: {} \n" \
          "Number of core pfams: {}".format(
              clusterID, 
              num_query_sublists + num_ref_sublists, 
              num_query_sublists, 
              num_ref_sublists,
              core_genome_indicator, 
              antismash_indicator,
              transporter_indicator, 
              number_core_pfams
              ), file=outfile)
    print("CDS products:", file=outfile)
    for cds_product in cds_in_cluster:
        print('   ', cds_product, file=outfile)
    print("", file=outfile)


def create_xlsx_writer(out_dir):
    detailed_results_path = os.path.join(out_dir, "results_detailed.xlsx")
    writer = pd.ExcelWriter(detailed_results_path, engine='xlsxwriter')
    workbook  = writer.book
    return writer, workbook


def get_results(cluster_list, database_path, core_genome_format, 
                antismash_format):

    conn = sqlite3.connect(database_path)

    summary_results = list()
    detailed_results = list()

    for cluster in cluster_list:
        
        clusterID = cluster[0]
        core_genome_indicator = cluster[1]
        antismash_indicator = cluster[2]
        transporter_indicator = cluster[3]
        number_core_pfams = cluster[4]

        # Get all the sublists in the cluster
        sublists = _get_sublists(conn, clusterID) 
        
        # For summary results
        summary_results.append([
            clusterID,
            core_genome_indicator,
            antismash_indicator,
            transporter_indicator,
            number_core_pfams,
            _num_sublist_of_host_type(sublists, "query"), 
            _num_sublist_of_host_type(sublists, "ref"),
            _labeled_summary_cds_products(sublists)
            ])

        # for detailed results
        for sublist in sublists: 
            detailed_results.append([
                clusterID,
                len(sublists), 
                core_genome_indicator, 
                antismash_indicator,
                transporter_indicator, 
                number_core_pfams,
                sublist.id,
                sublist.host_type, 
                sublist.organism, 
                sublist.file, 
                sublist.seq_acc, 
                sublist.description, 
                sublist.start, 
                sublist.end, 
                _labeled_cds_products_set(sublist),  
                _labeled_cds_products(sublist, core_genome_format, 
                                      antismash_format), 
                _labeled_updownstream(sublist, core_genome_format, 
                                      "upstream"),
                _labeled_updownstream(sublist, core_genome_format, 
                                      "downstream"),
                ])
    
    conn.close()
    
    return summary_results, detailed_results


def write_detailed_results(detailed_results, writer, workbook):
    
    # Write detailed results to a tsv file
    cols = [
        'clusterID', 
        'number of sublists', 
        'core genome indicator', 
        'antismash indicator',
        'transporter indicator', 
        'number of core pfams', 
        'sublist id',
        'host type', 
        'organism', 
        'file', 
        'sequence accession', 
        'sequence description', 
        'start on sequence', 
        'end on sequence', 
        'CDS products set', 
        'CDS products sequence', 
        'upstream CDS products', 
        'downstream CDS products',
        ]
    detailed_results = pd.DataFrame(detailed_results, columns = cols)
    detailed_results.to_excel(writer, sheet_name='Sheet1', index=False)
    worksheet = writer.sheets['Sheet1']
    
    # Add autofilter function
    worksheet.autofilter(0, 0, len(detailed_results), 
                         len(detailed_results.columns) - 1)
    
    # set column width 
    worksheet.set_column(first_col=cols.index('sequence description'), 
                         last_col =cols.index('sequence description'), 
                         width=40)
    worksheet.set_column(first_col=cols.index('CDS products set'), 
                         last_col =cols.index('downstream CDS products'), 
                         width=60)
    
    # bg color for ref hosts
    bg_color = "#b2d4b9"
    border_color = "#C0C0C0"
    ref_row_format = workbook.add_format(
        {'border': 1, 'bg_color': bg_color, 'border_color': border_color})
    
    for row_idx in range(0, len(detailed_results)):
        row_num = row_idx + 1, # worksheet is 1 index based
        # fill cells of ref host rows
        if detailed_results.iloc[row_idx, 6] == "ref":
            worksheet.set_row(
                row=row_num,
                height=None, 
                cell_format=ref_row_format)
            # format for cds cells
            cds_cells_format = workbook.add_format(
                {
                    'border': 1, 
                    'text_wrap': True, 
                    'border_color': border_color, 
                    'bg_color': bg_color,
                 }
                )
        else:
            # format for cds cells
            cds_cells_format = workbook.add_format(
                {
                    'border': 1, 
                    'text_wrap': True, 
                    'border_color': border_color}
                )

        # write the CDS as rich strings $$$
        cds_indeces = range(
            cols.index('CDS products sequence'),
            cols.index('downstream CDS products')+1
            )
        for col_idx in cds_indeces:
            cds = detailed_results.iloc[row_idx, col_idx]
            # write_rich_string() needs at least 3 fragments/formats
            if len(cds) >= 3: 
                worksheet.write_rich_string(row_idx+1, col_idx, *cds, 
                                            cds_cells_format)
            # no cds
            elif len(cds) == 0:
                worksheet.write_blank(row_idx+1, col_idx, "", 
                                      cds_cells_format)
            # only one cds without underline
            elif len(cds) == 1:
                worksheet.write_string(row_idx+1, col_idx, *cds, 
                                       cds_cells_format)
            # only one cds with underline
            elif len(cds) == 2:
                cds_cells_format.set_underline()
                worksheet.write_string(row_idx+1, col_idx, cds[1], 
                                       cds_cells_format)

    workbook.close()


def write_summary_file(summary_results, database_path, core_genome_cutoff, 
                        antismash_cutoff, transporter_cutoff, out_dir):

    # Write summary results to txt file
    summary_results_path = os.path.join(out_dir, "results_summary.txt")
    with open(summary_results_path, "w") as outfile_summary:
        _print_summary_header(
            database_path, 
            core_genome_cutoff, 
            antismash_cutoff, 
            transporter_cutoff, 
            outfile_summary
            )
        for entry in summary_results:
            _print_cluster_to_summary_results(entry, outfile_summary)
