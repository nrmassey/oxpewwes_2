#!/usr/bin/env python
###############################################################################
# Program : convert_ascii_pop_to_netcdf.py
# Author  : Neil Massey
# Date	  : 29/10/13
# Purpose : Convert the ascii population files to netCDF:
#
# Center for International Earth Science Information Network (CIESIN)/Columbia 
# University, and Centro Internacional de Agricultura Tropical (CIAT). 2005. 
# Gridded Population of the World, Version 3 (GPWv3): Subnational Administrative 
# Boundaries. Palisades, NY: NASA Socioeconomic Data and Applications Center 
# (SEDAC). http://sedac.ciesin.columbia.edu/data/set/gpw-v3-subnational-admin-boundaries. 
# Accessed 28 OCT 2013.
###############################################################################

import sys, os, getopt
from scipy.io.netcdf import *
import numpy

###############################################################################

def read_input_file(in_file):
	fh = open(in_file)
	lines = fh.readlines()
	# get the metadata from the file
	n_cols = int(lines[0].split()[1])
	n_rows = int(lines[1].split()[1])
	xll = float(lines[2].split()[1])
	yll = float(lines[3].split()[1])
	deg_per_cell = float(lines[4].split()[1])
	mv = float(lines[5].split()[1])
	# create the arrays
	lon = numpy.zeros((n_cols,), 'f')
	lat = numpy.zeros((n_rows,), 'f')	
	# put the lat / lon values in the array
	for c in range(0, n_cols):
		lon[c] = xll + c*deg_per_cell
	for r in range(0, n_rows):
		lat[r] = yll + r*deg_per_cell
	array = numpy.zeros([n_rows, n_cols], 'f')
	for r in range(0, n_rows):
		c_line = lines[r+6].split()
		for c in range(0, n_cols):
			array[n_rows-1-r,c] = float(c_line[c])
	fh.close()
	return array, lon, lat, mv

###############################################################################

def save_file(out_file, array, lon, lat, mv):
	# create output file
	out_file = netcdf_file(out_file, "w")
	# create the dimensions
	lon_dim  = out_file.createDimension("lon", lon.size)
	lat_dim  = out_file.createDimension("lat", lat.size)
	lon_var  = out_file.createVariable("lon", lon.dtype, ("lon",))
	lat_var  = out_file.createVariable("lat", lat.dtype, ("lat",))
	lon_var.__setattr__("standard_name", "longitude")
	lon_var.__setattr__("axis", "X")
	lon_var.__setattr__("units", "degrees_east")
	lat_var.__setattr__("standard_name", "latitude")
	lat_var.__setattr__("units", "degrees_north")
	# assign the data
	lon_var[:] = lon
	lat_var[:] = lat
	# create the variable
	data_var = out_file.createVariable("population", array.dtype, ("lat", "lon"))
	data_var[:,:] = array[:,:]
	# add the attributes to the variable
	data_var.__setattr__("missing_value", mv)
	data_var.__setattr__("units", "persons per km^2")
	data_var.__setattr__("origin", "http://sedac.ciesin.columbia.edu/data/set/gpw-v3-subnational-admin-boundaries")
	data_var.__setattr__("accessed", "2013-10-28")
	data_var.__setattr__("author", "Neil Massey, University of Oxford")
	out_file.close()

###############################################################################

if __name__ == "__main__":
	in_file  = ""
	out_file = ""
	
	opts, args = getopt.getopt(sys.argv[1:], 'i:o:', 
							   ['in_name=', 'out_name=',])

	for opt, val in opts:
		if opt in ['--in_name', '-i']:
			in_file = val
		if opt in ['--out_name', '-o']:
			out_file = val

	# read the file into a numpy array
	array, lon, lat, mv = read_input_file(in_file)
	save_file(out_file, array, lon, lat, mv)