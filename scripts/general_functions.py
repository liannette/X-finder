import sqlite3
import os.path
import time

def get_base_dir():
    '''
    '''
    dirname = os.path.dirname
    BASE_DIR = dirname(dirname(os.path.abspath(__file__)))
    return BASE_DIR

def connect_to_db(db_filename):
    '''
    Connects to a database in BASE_DIR/data/database
    If the database does not exist, it will be created
    :param db_file: database file
    :return: Connection object or None
    '''
    database = os.path.join(get_base_dir(), "data", "database", db_filename)
    conn = sqlite3.connect(database, timeout=20)
    return conn



def timepoint(comment, timelist):
    now = time.time()
    timelist.append((comment, now))


def showtimepoints(timelist):
    m = max([ len(timelist[i][0]) for i in range(1, len(timelist)) ])
    formatstring = "{:" + str(m+1) +"} {:.02f} seconds"
    for i in range(1, len(timelist)):
        print(formatstring.format(timelist[i][0]+':', timelist[i][1]-timelist[i-1][1]))
    print(formatstring.format('Total:', timelist[-1][1]-timelist[0][1]))