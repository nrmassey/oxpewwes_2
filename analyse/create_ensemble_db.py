from create_db_common import *
import sqlite3
import os

##############################################################################

if __name__ == "__main__":
    DB_NAME = "oxpewwes_ensemble.db" 
    # create the database connection
    create_db = True 
    if os.path.exists(DB_NAME):
        create_db = False

    conn = sqlite3.connect(DB_NAME)

    if create_db:
        create_db_schema(conn)

    # root path to OXPEWWES_2 data
    root_dir = "/group_workspaces/jasmin2/cpdn_rapidwatch/OXPEWWES_2/"

    # get the list of ensemble members over the years
    for yr in range(1989, 2011):
        rel_dir = "ensemble/events/"+str(yr)+"_"+str(yr+1)+"/"
        ensemble_list = os.listdir(root_dir+rel_dir)
        # loop over the ensemble members - could be up to 800 per year
        for e in ensemble_list:
            # get the list of files for each ensemble member - should be two
            file_list = os.listdir(root_dir+rel_dir+e)
            for f in file_list:
                # get the relative path to the file
                rel_path = rel_dir + e + "/" + f
                # get the umid from the file / ensemble member
                umid = e.split("_")[2]
                # process the file to ingest into database
                process_file_to_db(root_dir, rel_path, umid, conn)

    conn.close()
