###############################################################################
# Program : read_write_event.py
# Author  : Neil Massey
# Date    : 03/04/17
# Purpose : Functions to read in an event and write out an event
#           
###############################################################################

from netCDF4 import Dataset
import numpy


class Event:
    """Class containing event information from reading or creating"""

    def __init__(self):
        pass


def get_all_attributes(nc_obj):
    """Get all the attributes : values pairs as a dictionary"""
#    attrs = {}
#    for attr in nc_obj.ncattrs():
#        attrs[attr] = nc_obj.getncattr(attr)
    return nc_obj.__dict__


def put_all_attributes(nc_obj, attributes):
    """Write all attributes to a netCDF object (file or variable)"""
    nc_obj.setncatts(attributes)


def read_event(input_fname):
    """Read an event from the file and return an Event object.
       This is a more fully-featured version as it also reads in the metadata, ready to be written out again."""
    # create the event first
    event = Event()

    # open the netCDF file first
    nc_fh = Dataset(input_fname, 'r', format="NETCDF4")

    # get the number of events
    event.n_evts = len(nc_fh.dimensions["event"])
    
    # get the grid latitude and longitude and rotated pole info
    event.grid_latitude = nc_fh.variables["latitude"][:]
    event.grid_latitude_attributes = get_all_attributes(nc_fh.variables["latitude"])
    event.grid_longitude = nc_fh.variables["longitude"][:]
    event.grid_longitude_attributes = get_all_attributes(nc_fh.variables["longitude"])
    event.rotated_pole_attributes = get_all_attributes(nc_fh.variables["rotated_pole"])

    # get the track details
    event.track_time_attributes = get_all_attributes(nc_fh.variables["track_time"])
    event.track_latitude_attributes = get_all_attributes(nc_fh.variables["track_latitude"])
    event.track_longitude_attributes = get_all_attributes(nc_fh.variables["track_longitude"])
    # each track has a variable length - store as a list of numpy arrays
    event.track_time = []
    event.track_latitude = []
    event.track_longitude = []

    # store each track time, latitude, longitude tuple
    for e in range(0, event.n_evts):
        event.track_time.append(nc_fh.variables["track_time"][e])
        event.track_latitude.append(nc_fh.variables["track_latitude"][e])
        event.track_longitude.append(nc_fh.variables["track_longitude"][e])

    # now read in the variables and their attributes
    event.mslp = nc_fh.variables["mslp"][:]
    event.mslp_attributes = get_all_attributes(nc_fh.variables["mslp"])

    event.wind_max = nc_fh.variables["wind_max"][:]
    event.wind_max_attributes = get_all_attributes(nc_fh.variables["wind_max"])

    event.wind_gust = nc_fh.variables["wind_gust"][:]
    event.wind_gust_attributes = get_all_attributes(nc_fh.variables["wind_gust"])

    event.loss = nc_fh.variables["loss"][:]
    event.loss_attributes = get_all_attributes(nc_fh.variables["loss"])

    event.precip = nc_fh.variables["precip"][:]
    event.precip_attributes = get_all_attributes(nc_fh.variables["precip"])

    # get the global attributes
    event.attributes = get_all_attributes(nc_fh)
    nc_fh.close()

    return event
    

def write_event(output_fname, event):
    """Write an event out to a netCDF4 file"""
    # create the file
    out_fh = Dataset(output_fname, "w")    

    # write out the global attributes
    put_all_attributes(out_fh, event.attributes)

    # write out the grid latitudes and longitudes
    out_fh.createDimension("latitude", event.grid_latitude.shape[0])
    out_fh.createDimension("longitude", event.grid_longitude.shape[0])
    lat_var = out_fh.createVariable("latitude", 'f4', ("latitude"), shuffle=True, zlib=True, complevel=5)
    lon_var = out_fh.createVariable("longitude", 'f4', ("longitude"), shuffle=True, zlib=True, complevel=5)
    put_all_attributes(lat_var, event.grid_latitude_attributes)
    put_all_attributes(lon_var, event.grid_longitude_attributes)
    # create the events dimension
    out_fh.createDimension("event", event.n_evts)

    # write the longitude and latitude values
    lat_var[:] = event.grid_latitude[:]
    lon_var[:] = event.grid_longitude[:]

    # create the rotated grid variable
    rot_var = out_fh.createVariable("rotated_pole", 'S1', ())
    put_all_attributes(rot_var, event.rotated_pole_attributes)

    # create the variable length for the tracks
    track_len_type = out_fh.createVLType('f4', "trackLen")

    # add the track latitude, longitude and time variables
    track_time_var = out_fh.createVariable("track_time", track_len_type, ("event"), shuffle=True, zlib=True, complevel=5)
    track_latitude_var = out_fh.createVariable("track_latitude", track_len_type, ("event"), shuffle=True, zlib=True, complevel=5)
    track_longitude_var = out_fh.createVariable("track_longitude", track_len_type, ("event"), shuffle=True, zlib=True, complevel=5)
    put_all_attributes(track_time_var, event.track_time_attributes)
    put_all_attributes(track_latitude_var, event.track_latitude_attributes)
    put_all_attributes(track_longitude_var, event.track_longitude_attributes)

    # add the field variables
    # mslp
    mslp_var = out_fh.createVariable("mslp", 'i1', ("event", "latitude", "longitude"), shuffle=True, zlib=True, complevel=5)
    put_all_attributes(mslp_var, event.mslp_attributes)
    # wind_max
    wind_max_var = out_fh.createVariable("wind_max", 'i1', ("event", "latitude", "longitude"), shuffle=True, zlib=True, complevel=5)
    put_all_attributes(wind_max_var, event.wind_max_attributes)
    # wind_gust
    wind_gust_var = out_fh.createVariable("wind_gust", 'i1', ("event", "latitude", "longitude"), shuffle=True, zlib=True, complevel=5)
    put_all_attributes(wind_gust_var, event.wind_gust_attributes)
    # loss
    loss_var = out_fh.createVariable("loss", 'i1', ("event", "latitude", "longitude"), shuffle=True, zlib=True, complevel=5)
    put_all_attributes(loss_var, event.loss_attributes)
    # precip
    precip_var = out_fh.createVariable("precip", 'i1', ("event", "latitude", "longitude"), shuffle=True, zlib=True, complevel=5)
    put_all_attributes(precip_var, event.precip_attributes)

    # write the fields out, require packing first
    # mslp
    mslp_var[:] = event.mslp
    # wind_max
    wind_max_var[:] = event.wind_max
    # wind_gust
    wind_gust_var[:] = event.wind_gust
    # loss
    loss_var[:] = event.loss
    # precip
    precip_var[:] = event.precip

    # write the tracks out
    for e in range(0, event.n_evts):
        track_time_var[e] = event.track_time[e]
        track_latitude_var[e] = event.track_latitude[e]
        track_longitude_var[e] = event.track_longitude[e]

    out_fh.close()

if __name__ == "__main__":
    # test - read in an event and write it straight back out
    test_fname = "/group_workspaces/jasmin2/cpdn_rapidwatch/OXPEWWES_2/calibration/events/1989_1990/oxfaga_1989-12-01T01-00-00_1990-04-01T19-00-00.nc"
    output_fname = "./oxfaga_1989-12-01T01-00-00_1990-04-01T19-00-00.nc"

    event = read_event(test_fname)
    write_event(output_fname, event)
