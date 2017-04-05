###############################################################################
# Program : create_db_common.py
# Author  : Neil Massey
# Date    : 03/02/17
# Purpose : Common functions for creating the databases of the wind-storm event
#           footprints
###############################################################################

import numpy
from netCDF4 import Dataset
import sqlite3
import sys
sys.path.append("../plot")
from load_data import load_data
from get_date import get_date

###############################################################################

def create_db_schema(conn):
    """Create the initial database schema.  This is the same for both the
       calibration and ensemble event set."""
    # get a cursor to the database
    c = conn.cursor()

    # create the schema / table (see the file db_schema.txt for info)
    db_str = '''CREATE TABLE event
                 (umid TEXT, file_path TEXT, event_index INTEGER, 
                  start_date_year INTEGER, start_date_month INTEGER, start_date_day INTEGER, start_date_hour INTEGER,
                  end_date_year INTEGER, end_date_month INTEGER, end_date_day INTEGER, end_date_hour INTEGER,
                  persistence INTEGER,
                  wind_max REAL, wind_max_land REAL,
                  gust_max REAL, gust_max_land REAL,
                  precip_max REAL, precip_max_land REAL, precip_sum_land REAL,
                  loss_max REAL, loss_sum REAL,
                  mslp_min REAL)'''
    c.execute(db_str)
    conn.commit()

###############################################################################

def process_file_to_db(root_dir, relative_path, umid, conn):
    """Process data from the netcdf file containing the events
       inputs: root_dir - the first part of the file name, where the event set is kept
               relative_path - the path relative to the root dir where the event is
               umid - the UM id of the run that generated the event set
               conn - the database connection"""

    # get a cursor to the database
    c = conn.cursor()

    # load the land sea mask
    lsm_fh = Dataset("../wah_eu_lsm_0.44.nc")
    lsm_var = lsm_fh.variables["lsm"]
    lsm = lsm_var[:].squeeze()
    lsm_fh.close()

    # open the netCDF file first
    nc_fh = Dataset(root_dir+relative_path, 'r', format="NETCDF4")

    # get the number of events
    n_evts = len(nc_fh.dimensions["event"])

    # get the variables required
    time_var = nc_fh.variables["track_time"]
    time_units = time_var.units
    calendar_type = time_var.calendar
    wind_var = nc_fh.variables["wind_max"]
    gust_var = nc_fh.variables["wind_gust"]
    mslp_var = nc_fh.variables["mslp"]
    loss_var = nc_fh.variables["loss"]
    precip_var = nc_fh.variables["precip"]

    # loop over each event and read in the data
    for e in range(0, n_evts):
        # get the start and end time and convert to years-months-days:hrs
        time_start = time_var[e][0]
        time_end   = time_var[e][-1]

        time_start_year, time_start_month, time_start_day, time_start_hour = get_date(time_start, time_units, calendar_type)
        time_end_year, time_end_month, time_end_day, time_end_hour = get_date(time_end, time_units, calendar_type)

        # get the persistence
        persistence = int((time_end - time_start) * 24 + 0.5)
        
        # get the field data
        wind = wind_var[e][:]
        gust = gust_var[e][:]
        mslp = mslp_var[e][:]
        loss = loss_var[e][:]
        precip = precip_var[e][:] * 60.0 * 60.0 # (convert to mm/hr)
        
        # get the max and min values
        wind_max = numpy.max(wind)
        gust_max = numpy.max(gust)
        mslp_min = numpy.min(mslp)
        loss_max = numpy.max(loss)
        loss_sum = numpy.sum(loss)
        precip_max = numpy.max(precip)

        # get the max and min values over land
        wind_max_land = numpy.max(wind * lsm)
        gust_max_land = numpy.max(gust * lsm)
        precip_max_land = numpy.max(precip * lsm)
        precip_sum_land = numpy.sum(precip * lsm)

        if wind_max is numpy.ma.masked or gust_max is numpy.ma.masked or mslp_min is numpy.ma.masked or \
           loss_max is numpy.ma.masked or loss_sum is numpy.ma.masked or precip_max is numpy.ma.masked:
            continue

        # insert into the database
        insert_string = "INSERT INTO event VALUES (\'" + umid + "\', "+\
                        "'" + relative_path + "', " + str(e) + ", " +\
                        str(time_start_year) + ", " + str(time_start_month) + ", " + str(time_start_day) + ", " + str(time_start_hour) + ", "+\
                        str(time_end_year) + ", " + str(time_end_month) + ", " + str(time_end_day) + ", " + str(time_end_hour) + ", " +\
                        str(persistence) + ", " +\
                        str(wind_max) + ", " + str(wind_max_land) + ", " +\
                        str(gust_max) + ", " + str(gust_max_land) + ", " +\
                        str(precip_max) + ", " + str(precip_max_land) + ", " + str (precip_sum_land) + ", " +\
                        str(loss_max) + ", " + str(loss_sum) + ", " +\
                        str(mslp_min) + ")"
        c.execute(insert_string)
    conn.commit()

    nc_fh.close()

