from create_db_common import *
import sqlite3
import os

##############################################################################

if __name__ == "__main__":
    DB_NAME = "oxpewwes_calibration.db"
    # create the database connection
    create_db = True
    if os.path.exists(DB_NAME):
        create_db = False
    
    conn = sqlite3.connect(DB_NAME)

    if create_db:
        create_db_schema(conn)

    # root path to OXPEWWES_2 data
    root_dir = "/group_workspaces/jasmin2/cpdn_rapidwatch/OXPEWWES_2/"

    # umid for calibration data is "oxfa"
    umid = "oxfa"

    # get the list of files over the years
    for yr in range(1989, 2009):
        rel_dir = "calibration/events/"+str(yr)+"_"+str(yr+1)+"/"
        file_list = os.listdir(root_dir+rel_dir)
        # loop over the files (should only be one in calibration db)
        for f in file_list:
            # get the relative path to the file
            rel_path = rel_dir + f
            # process the file to ingest into database
            process_file_to_db(root_dir, rel_path, umid, conn)
    
    conn.close()
