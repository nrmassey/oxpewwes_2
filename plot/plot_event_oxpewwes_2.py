#!/usr/bin/env python

###############################################################################
# Program : plot_event_oxpewwes_2.py
# Author  : Neil Massey
# Date    : 29/08/16
# Purpose : Program to plot an oxpewwes 2 event footprint
#           Generates a four panel plot:
#                MSLP | Precip
#               ------+-------
#                Gust | Loss
#           Each plot also has the storm track on it
###############################################################################

import matplotlib.pyplot as plt
import matplotlib.colors as col
import matplotlib.cm as cm
import matplotlib
import sys, os, getopt
import cartopy.crs as ccrs
import numpy
from scipy.io.netcdf import *

###############################################################################

def create_wind_color_map(levels):
    # replicate XWS colours for wind speed
    #
    cmap = ["#ffffff", "#fff4e8", "#ffe1c1", "#ffaa4e", "#ff6d00",
            "#d33100", "#893030", "#652828", "#502020", "#401818", 
            "#301010", "#100000"]
    ccmap, norm = col.from_levels_and_colors(levels, cmap[:len(levels)-1], 'neither')
    return ccmap, norm

###############################################################################

def create_precip_color_map(levels):
    # colour map for precip amount - once scaled to mm/day
    cmap = ["#dddddd", "#0000ff", "#00ffff", "#ffff00", 
            "#ff0000", "#ff00ff"]
    ccmap, norm = col.from_levels_and_colors(levels, cmap, 'neither')
    return ccmap, norm
    
###############################################################################

def create_mslp_color_map(levels):
    # colour map for mslp amount
    cmap = ["#003366", "#0055bb", "#0088ff", "#00ffff", 
            "#ffff00", "#ff8800", "#ff0000", "#880000"]
    ccmap, norm = col.from_levels_and_colors(levels, cmap, 'neither')
    return ccmap, norm

###############################################################################

def create_loss_color_map(levels):
    # colour map fo loss amount
    cmap = ["#dddddd", "#a5ffff", "#00ff00", "#ffff00", "#ffbf00", "#ff0000"]
    ccmap, norm = col.from_levels_and_colors(levels, cmap, 'neither')
    return ccmap, norm

###############################################################################

def plot_track(sp0, track_lons, track_lats):
    # plot the track, first plot markers at each 6h timestep
    for i in range(0, track_lons.shape[0]):
        sp0.plot(track_lons[i], track_lats[i], 'ko', ms=3, 
                 transform=ccrs.PlateCarree())
    # plot the line between the points
    for i in range(0, track_lons.shape[0]-1):
        lon_0 = track_lons[i]
        lon_1 = track_lons[i+1]
        
        if lon_0 > 180.0 and lon_1 < 180.0:
            lon_0 -= 360.0
        if lon_1 > 180.0 and lon_0 < 180.0:
            lon_1 -= 360.0
        
        sp0.plot((lon_0, lon_1), (track_lats[i], track_lats[i+1]),
                 'k', lw=2.0, transform=ccrs.PlateCarree())

###############################################################################

def load_scaled_data(nc_fh, var_name):
    # load data that has been packed to a byte or short and scaled
    var = nc_fh.variables[var_name]
    # get the scale factor and offset
    sf = var._attributes["scale_factor"]
    off = var._attributes["add_offset"]
    # get the missing value
    mv = var._attributes["_FillValue"]
    # get the data and mask
    D = numpy.ma.masked_equal(var[:], mv)
    # scale and return
    SD = D * sf + off
    return SD

###############################################################################

if __name__ == "__main__":
    
    font = {'family' : 'sans-serif',
            'weight' : 'medium',
            'size'   : 12}

    matplotlib.rc('font', **font)

    input_file = ""
    output_file = ""
    # get the input and output filenames
    opts, args = getopt.getopt(sys.argv[1:], 'i:o:l:', 
                               ['input=', 'output=', 'min_loss='])

    for opt, val in opts:
        if opt in ['--input', '-i']:
            input_file = val
        if opt in ['--output', '-o']:
            output_file = val
        if opt in ['--min_loss', '-l']:
            min_loss = float(val)

    if input_file == "" or output_file == "":
        print " Usage : plot_event_oxpewwes_2 -i input -o output "
        sys.exit()

    # load the input file
    fh = netcdf_file(input_file, 'r', mmap=False)
    # get the pole latitude and longitude from the "rotated_pole" variable
    rot_grid_var = fh.variables["rotated_pole"]
    rot_lonp = float(rot_grid_var._attributes["grid_north_pole_longitude"])
    rot_latp = float(rot_grid_var._attributes["grid_north_pole_latitude"])
    
    # get the data
    mslp_data = load_scaled_data(fh, "mslp") * 0.01
    gust_data = load_scaled_data(fh, "wind_gust")
    precip_data = load_scaled_data(fh, "precip") * 60.0 * 60.0  # convert to mm/hr
    loss_data = load_scaled_data(fh, "loss")
    lons = fh.variables["longitude"][:]
    lats = fh.variables["latitude"][:]
    
    # get the loss amount
    loss_val = float(fh.variables["loss"]._attributes["sum"])
    if loss_val < min_loss:
        sys.exit()
    
    # get the track data
    track_lons = fh.variables["track_longitude"][:]
    track_lats = fh.variables["track_latitude"][:]
    
    # create the projection
    proj = ccrs.RotatedPole(pole_latitude=rot_latp, pole_longitude=rot_lonp)

    # create the subplots
    sps = []
    sps.append(plt.subplot(221, projection=proj))
    sps.append(plt.subplot(222, projection=proj))
    sps.append(plt.subplot(223, projection=proj))
    sps.append(plt.subplot(224, projection=proj))
    
    # create the color maps and levels
    wind_levels = [x for x in range(0,65,5)]
    precip_levels = [0, 0.1, 1, 2, 4, 6, 8]
    mslp_levels = [x for x in range(910,1090,20)]
    loss_levels = [0, 5, 10, 20, 30, 40, 50]
    
    wind_cmap, wind_norm = create_wind_color_map(wind_levels)
    precip_cmap, precip_norm = create_precip_color_map(precip_levels)
    mslp_cmap, mslp_norm = create_mslp_color_map(mslp_levels)
    loss_cmap, loss_norm = create_loss_color_map(loss_levels)
    
    # mslp
    mslp_map = sps[0].pcolormesh(lons, lats, mslp_data,
                                 cmap=mslp_cmap, norm=mslp_norm,
                                 vmax=mslp_levels[-1], vmin=mslp_levels[0])
    sps[0].set_title("MSLP hPa")

    # draw colorbar
    mslp_bar = plt.gcf().colorbar(mslp_map, ax = sps[0], orientation="horizontal",
                                  fraction=0.05, pad=0.04)
#    mslp_bar.set_label("MSLP hPa")
    
    # precip
    precip_map = sps[1].pcolormesh(lons, lats, precip_data, 
                                   cmap=precip_cmap, norm=precip_norm,
                                   vmax=precip_levels[-1], vmin=precip_levels[0])
    sps[1].set_title("Precip. mm hr$^{-1}$")
    precip_bar = plt.gcf().colorbar(precip_map, ax = sps[1], orientation="horizontal",
                                    fraction=0.05, pad=0.04)
#    precip_bar.set_label("Precip. mm hr$^{-1}$")

    # wind
    wind_map = sps[2].pcolormesh(lons, lats, gust_data,
                                 cmap=wind_cmap, norm=wind_norm,
                                 vmax=wind_levels[-1], vmin=wind_levels[0])
    sps[2].set_title("3s wind gust (m s$^{-1}$)")
    wind_bar = plt.gcf().colorbar(wind_map, ax = sps[2], orientation="horizontal",
                                  fraction=0.05, pad=0.04)
#    wind_bar.set_label("3s wind gust (m s$^{-1}$)")

    loss_map = sps[3].pcolormesh(lons, lats, loss_data,
                                 cmap=loss_cmap, norm=loss_norm,
                                 vmax=loss_levels[-1], vmin=loss_levels[0])
    sps[3].set_title("Est. losses")
    loss_bar = plt.gcf().colorbar(loss_map, ax = sps[3], orientation="horizontal",
                                  fraction=0.05, pad=0.04)
#    loss_bar.set_label("Est. losses")
    
    # draw the continents etc. on the subplots
    for sp in sps:
        # plot the track on each subplot
        plot_track(sp, track_lons, track_lats)
        sp.coastlines(resolution='50m', lw=0.5, zorder=3)
        sp.gridlines()
        sp.get_axes().set_extent([-15.5, 39.0, 28.5, 72.0])
#        sp.stock_img()
        sp.set_aspect(1.0)

    # save the figure
    plt.gcf().set_size_inches(9,9)
    plt.tight_layout()
    plt.savefig(output_file)
    
    # close the file
    fh.close()