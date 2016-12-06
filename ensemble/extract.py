#!/usr/bin/env python

###############################################################################
# Program : wah_extract.py
# Author  : Neil Massey
# Date    : 04/06/14
# Purpose : download w@h output files from a list of urls and extract the 
#           requested fields into separate netCDF files
# Requires: scipy.io.netcdf
###############################################################################

import sys, getopt, os, ast
import urllib, tempfile, zipfile, shutil
import numpy, numpy.ma
import math
from scipy.io.netcdf import *

###############################################################################

def strip(s):
    return s.strip()

###############################################################################

def read_urls(urls_file):
    fh = open(urls_file, 'r')
    urls = fh.readlines()
    fh.close()
    urls = map(strip, urls)
    return urls

###############################################################################

def get_umid_boinc_year(url):
    # get the umid, boinc id and year from url
    # get the last file
    lf = url.split("/")[-1][:-4]
    boinc = url.split("/")[-2]
    umid = lf.split("_")[2]
    year = lf.split("_")[3]
    return umid, boinc, year

###############################################################################

def get_output_field_name(field):
    fname = field[1]
    return fname

###############################################################################

def get_output_dir(output_dir, year, boinc):
    return output_dir + "/" + str(year) + "/" + boinc

###############################################################################

def make_directories(output_dir, year, boinc, field_list, n_valid):
    # make the output top directory
    c_path = output_dir
    download=False
    if not os.path.exists(c_path):
        os.mkdir(c_path)
    # make the year directory
    c_path += "/" + year
    if not os.path.exists(c_path):
        os.mkdir(c_path)
    # make the boinc / umid directory
    c_path += "/" + boinc
    if not os.path.exists(c_path):
        os.mkdir(c_path)
    # finally - make the field directories
    for field in field_list:
        # make the pattern directory (i.e. pattern for global / regional / regional daily etc.)
        f_path = c_path + "/" + field[0]
        if not os.path.exists(f_path):
            os.mkdir(f_path)
        # if there is not subsetting or processing of the field then the name is 
        # just the field name
        f_path = f_path + "/" + get_output_field_name(field)
            
        if not os.path.exists(f_path):
            os.mkdir(f_path)
        if len(os.listdir(f_path)) != n_valid:
            download=True

    return c_path, download

###############################################################################

def get_missing_value(attrs):
    if "missing_value" in attrs.keys():
        mv = attrs["missing_value"]
    if "_FillValue" in attrs.keys():
        mv = attrs["_FillValue"]
    return mv

###############################################################################

def extract_netcdf_var(nc_in_file, nc_out_file, field):
    nc_in_var = nc_in_file.variables[field[1]]
    in_dimensions = []

    # now copy the dimensions from input netcdf
    for d in nc_in_var.dimensions:
        # get the input dimension and the data
        dim_in_var = nc_in_file.variables[d]
        dim_in_data = dim_in_var[:]
        in_dimensions.append([d, dim_in_data])
    out_dims = in_dimensions

    # rename the output dimensions
    out_dim_names = []
    
    # create the variable
    for d in out_dims:
        # create the output dimension and variable
        if "time" in d[0]:
            out_name = "time"
        elif "latitude" in d[0]:
            out_name = "latitude"
        elif "longitude" in d[0]:
            out_name = "longitude"
        elif "z" in d[0]:
            out_name = "z"
        out_dim_names.append(out_name)
        nc_out_file.createDimension(out_name, d[1].shape[0])
        dim_out_var = nc_out_file.createVariable(out_name, d[1].dtype, (out_name,))
        # assign the output variable data and attributes from the input
        if d[0] in nc_in_file.variables.keys():
            dim_in_var = nc_in_file.variables[d[0]]
            dim_out_var._attributes = dim_in_var._attributes
        dim_out_var[:] = d[1][:]

    # copy the data
    var_out_data = nc_in_var[:,:,:,:]

    nc_out_var = nc_out_file.createVariable(field[1], var_out_data.dtype, out_dim_names)
    # assign the attributes
    nc_out_var._attributes = nc_in_var._attributes
    # assign the data
    nc_out_var[:] = var_out_data
    # check for rotated pole and copy variable if it exists
    if "grid_mapping" in nc_out_var._attributes and len(out_dims) == 4:
        grid_map_name = nc_out_var._attributes["grid_mapping"]
        grid_map_var = nc_in_file.variables[grid_map_name]
        grid_map_out_var = nc_out_file.createVariable(grid_map_name, 'c', ())
        grid_map_out_var._attributes = grid_map_var._attributes
        # get the global longitude / global latitude vars
        coords = (nc_out_var._attributes["coordinates"]).split(" ");
        global_lon_var = nc_in_file.variables[coords[0]]
        global_lat_var = nc_in_file.variables[coords[1]]
        global_lon_data = global_lon_var[:,:]
        global_lat_data = global_lat_var[:,:]
        # create the global latitude / global longitude variables
        out_global_lon_var = nc_out_file.createVariable(coords[0], global_lon_data.dtype, (out_dim_names[2], out_dim_names[3]))
        out_global_lon_var[:] = global_lon_data
        out_global_lat_var = nc_out_file.createVariable(coords[1], global_lat_data.dtype, (out_dim_names[2], out_dim_names[3]))
        out_global_lat_var[:] = global_lat_data
        out_global_lon_var._attributes = global_lon_var._attributes
        out_global_lat_var._attributes = global_lat_var._attributes
        
###############################################################################

def extract_netcdf(zf, zf_file, base_path, field, temp_dir):
    # extract the file to a tempfile
    tmp_file_path = temp_dir + "/" + zf_file
    zf.extract(zf_file, temp_dir)
    # open as netCDF to a temporary file
    nc_file = netcdf_file(tmp_file_path)
    # create the output netCDF file
    o_field_name = get_output_field_name(field)
    out_name = base_path + "/" + field[0] + "/" + o_field_name + "/" + zf_file[:-3] + "_" + o_field_name + ".nc"
    # check whether it exists
    if os.path.exists(out_name):
        return
    # get the variable from the input
    if not field[1] in nc_file.variables.keys():
        print "Could not extract field: " + field[1] + " from file: " + zf_file
        return
    out_ncf = netcdf_file(out_name, "w")
    extract_netcdf_var(nc_file, out_ncf, field)
    out_ncf.close()
    nc_file.close()
    # delete the temp directory and its contents
    os.remove(tmp_file_path)

###############################################################################

def extract(url, field_list, output_dir, temp_dir, n_valid):
    # fetch the url to a temp (zip) file using urllib
   print url
   if True:
#   try:
        # get the umid, boinc id and year from the url
        umid, boinc, year = get_umid_boinc_year(url)
        # make directories to hold the output
        base_path, download = make_directories(output_dir, year, boinc, field_list, n_valid)
        # check whether this has already been processed
        if not download:
            return
        zf_fh = tempfile.NamedTemporaryFile(mode='w', delete=False)
        urlret = urllib.urlretrieve(url, zf_fh.name)
        # open the zip
        zf = zipfile.ZipFile(zf_fh.name)
        # list the zip contents
        zf_list = zf.namelist()
        # check to see which files match the pattern
        for field in field_list:
            pattern = field[0]
            for zf_file in zf_list:
                if pattern in zf_file:
                    extract_netcdf(zf, zf_file, base_path, field, temp_dir)
        os.remove(zf_fh.name)
   else:
#   except:
        print "Could not extract url: " + url

###############################################################################