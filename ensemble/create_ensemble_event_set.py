#! /usr/bin/env python

###############################################################################
# Program: create_ensemble_event_set.py
# Purpose: Build the OxPEWWES 2 ensemble event set
# Author : Neil Massey
# Date   : 12/11/2016
# Contact: neil.massey@ouce.ox.ac.uk
###############################################################################

from extract import *
import shutil
from cdo import *
import glob
import subprocess, os

###############################################################################

def read_batch_file(name, in_year):
    # Read in the batch file containing the list of ensemble members and
    # return a list of url prefixs
    fh = open(name)
    # read each line in from the file
    ls = fh.readlines()
    # ensemble members
    esm = []
    # split each line
    for l in ls:
        umid, boinc, year = get_umid_boinc_year(l)
        if in_year != year:
            continue
        s = l.split("/")
        # workunit name is in 6th place
        if len(s) > 6 and not s[6] in esm:
            esm.append(s[6])
    fh.close()
    return esm

###############################################################################

def get_url_prefix():
    p = "http://cpdn-upload2.oerc.ox.ac.uk/results/hadam3p_eu/batch48/"
    return p

###############################################################################

def get_urls_for_batch(b):
    # get the urls for Dec->Apr and Sep->Oct
    urls = []
    for i in range(1,6):
        url = get_url_prefix() + b + "/" + b + "_" + str(i) + ".zip"
        urls.append(url)
    for i in range(10,13):
        url = get_url_prefix() + b + "/" + b + "_" + str(i) + ".zip"
        urls.append(url)
    return urls

###############################################################################

def concatenate(u, field_list, output_dir):
    # get the actual output directory
    umid, boinc, year = get_umid_boinc_year(u)
    output_path = get_output_dir(output_dir, year, boinc)
    # months for first file
    mon_set_1 = ['dec', 'jan', 'feb', 'mar', 'apr']
    # months for second file
    mon_set_2 = ['sep', 'oct', 'nov']
    # cdo instance
    cdo = Cdo()
    # loop over the field_list
    for f in field_list:
        this_path = output_path + "/" + f[0] + "/" + f[1]
        # output file names
        output_fname_1 = this_path + "/" + umid + "_" + f[1] + "_DJFMA.nc"
        output_fname_2 = this_path + "/" + umid + "_" + f[1] + "_SON.nc"
        # only do this if the files don't already exist
        if os.path.exists(output_fname_1):
            continue
        mon_files_1 = []
        # build the first set of months
        for m in mon_set_1:
            fname = glob.glob(this_path + "/*" + m + "*.nc")
            if len(fname) != 0:
                mon_files_1.append(fname[0])
                # fix the record dimension using ncks
                cmd = ["ncks", "-O", "--mk_rec_dmn", "time", fname[0], "-o", fname[0]]
#                subprocess.call(cmd)
        # concat together using ncrcat
        cmd = ["ncrcat"]
        cmd.extend(mon_files_1)
        cmd.append(output_fname_1)
#        subprocess.call(cmd)
        
        # do the SON file          
        if os.path.exists(output_fname_2):
            continue
        mon_files_2 = []
        # build the second set of months
        for m in mon_set_2:
            fname = glob.glob(this_path + "/*" + m + "*.nc")
            if len(fname) != 0:
                mon_files_2.append(fname[0])
                # fix the record dimension using ncks
                cmd = ["ncks", "-O", "--mk_rec_dmn", "time", fname[0], "-o", fname[0]]
#                subprocess.call(cmd)
        # concat together using ncrcat
        cmd = ["ncrcat"]
        cmd.extend(mon_files_2)
        cmd.append(output_fname_2)
#        subprocess.call(cmd)
    return output_fname_1, output_fname_2, output_path

###############################################################################

def regrid(DJFMA_fname, SON_fname, LEV):

    # variables from the old scripts - field8 is MSLP
    N_PTS = 0
    PMODE = "0"
    GRID_FILE = "./grids/wah_mesh_EU_L"+str(LEV)
    NC_VAR = "field8"
    EXE_PATH=os.path.expanduser("~/Coding/tri_tracker/exe/regrid")

    # output file names
    DJFMA_regrid_fname = DJFMA_fname[:-3] + "_L" + str(LEV) + ".rgd"
    SON_regrid_fname = SON_fname[:-3] + "_L" + str(LEV) + ".rgd"
    
    # run the command - DJFMA then SON
    cmd = [EXE_PATH, "-i", DJFMA_fname, "-v", NC_VAR, "-m", GRID_FILE, "-o", DJFMA_regrid_fname, "-p", PMODE]
#    subprocess.call(cmd)
    cmd = [EXE_PATH, "-i", SON_fname, "-v", NC_VAR, "-m", GRID_FILE, "-o", SON_regrid_fname, "-p", PMODE]
#    subprocess.call(cmd)
    
    return DJFMA_regrid_fname, SON_regrid_fname

###############################################################################

def extrema(DJFMA_regrid_name, SON_regrid_name, LEV, EX_LEV, LS_LEV):

    PMODE = "0"
    GRID_FILE = "./grids/wah_mesh_EU_L"+str(LEV)
    EXE_PATH=os.path.expanduser("~/Coding/tri_tracker/exe/extrema")

    # build the filenames for the geopotential height
    GPH_DJFMA_fname = (DJFMA_regrid_name[:-7] + ".nc").replace("field8", "field1")
    GPH_SON_fname = (SON_regrid_name[:-7] + ".nc").replace("field8", "field1")
    
    # build the output filenames
    EX_DJFMA_fname = DJFMA_regrid_name[:-4] + "_E"+str(EX_LEV)+"_S"+str(LS_LEV)+".ex"
    EX_SON_fname = SON_regrid_name[:-4] + "_E"+str(EX_LEV)+"_S"+str(LS_LEV)+".ex"
    
    # build the arguments
    EX_DJFMA_args = "minima_largescale("+str(LS_LEV)+",100,-200)"
    GP_DJFMA_args = "geostrophic("+GPH_DJFMA_fname+",field1,0)"
    
    EX_SON_args = "minima_largescale("+str(LS_LEV)+",100,-200)"
    GP_SON_args = "geostrophic("+GPH_SON_fname+",field1,0)"
    
    # build the command - DJFMA
    cmd = [EXE_PATH, "-i", DJFMA_regrid_name, "-o", EX_DJFMA_fname, "-m", GRID_FILE, \
           "-l", str(EX_LEV), "-a", str(0), "--method", EX_DJFMA_args, \
           "--steering", GP_DJFMA_args]
#    subprocess.call(cmd)
    
    # SON
    cmd = [EXE_PATH, "-i", SON_regrid_name, "-o", EX_SON_fname, "-m", GRID_FILE, \
           "-l", str(EX_LEV), "-a", str(0), "--method", EX_SON_args, \
           "--steering", GP_SON_args]
#    subprocess.call(cmd)
    
    return EX_DJFMA_fname, EX_SON_fname

###############################################################################

def track(EX_DJFMA_fname, EX_SON_fname):

    EXE_PATH=os.path.expanduser("~/Coding/tri_tracker/exe/track")

    # build the output filenames
    TRK_DJFMA_fname = EX_DJFMA_fname[:-3]+".trk"
    TRK_SON_fname = EX_SON_fname[:-3]+".trk"
    
    # build the arguments
    TRK_ARGS = ["-r", "1000", "-f", "6", "-v", "10", "-s", "2"]
    cmd = [EXE_PATH, "-i", EX_DJFMA_fname, "-o", TRK_DJFMA_fname]
    cmd.extend(TRK_ARGS)
#    subprocess.call(cmd)
    
    cmd = [EXE_PATH, "-i", EX_SON_fname, "-o", TRK_SON_fname]
    cmd.extend(TRK_ARGS)
#    subprocess.call(cmd)

    return TRK_DJFMA_fname, TRK_SON_fname

###############################################################################

def events(DJFMA_track_name, SON_track_name, LEV, EX_LEV, LS_LEV, umid, boinc, year):
    # Build the event footprints from the tracks and input data
    EXE_PATH=os.path.expanduser("~/Coding/tri_tracker/exe/event_set")

    # netcdf variable names
    MSLP_VNAME = "field8"
    WIND_VNAME = "field50"
    PRECIP_VNAME = "field90"
    
    POP_VNAME = "population"
    LSM_VNAME = "lsm"
    
    # get the umid from the track file
    umid = DJFMA_track_name.split("/")[-1][0:4]
    output_path = "/".join(DJFMA_track_name.split("/")[0:-3])
    output_event_path = "/Volumes/SSD2/oxpewwes_2/ensemble/events/"
    
    # get the static input files for the loss / wind / population model
    POP_FILE=os.path.abspath("../population/euds00ag_wah_50km_final.nc")
    REMAP_FILE=os.path.abspath("../max_to_peak_wind.nc")
    LSM_FILE=os.path.abspath("../wah_eu_lsm_0.44.nc")
    
    # event output file
    # form the output_path
    # form the year directory
    event_path = output_event_path + str(year) + "_" + str(int(year)+1) 
    # add the boinc filename
    event_path += "/" + boinc + "/"
    if not os.path.exists(event_path):
        os.makedirs(event_path)
        
    EVENT_umid = event_path + umid
        
    # get the dynamic files - do the DJFMA files first
    MSLP_DJFMA = output_path + "/ga.pd/field8/" + umid + "_field8_DJFMA.nc"
    WIND_DJFMA = output_path + "/ga.pd/field50/" + umid + "_field50_DJFMA.nc"
    PRECIP_DJFMA = output_path + "/ga.pd/field90/" + umid + "_field90_DJFMA.nc"


    cmd = [EXE_PATH, "-m", MSLP_DJFMA, "-M", MSLP_VNAME,
                     "-w", WIND_DJFMA, "-W", WIND_VNAME,
                     "-p", PRECIP_DJFMA, "-P", PRECIP_VNAME,
                     "-q", POP_FILE, "-Q", POP_VNAME,
                     "-e", REMAP_FILE, "-l", LSM_FILE, "-L", LSM_VNAME,
                     "-t", "7", "-r", "1000.0", "-i",
                     DJFMA_track_name, "-o", EVENT_umid]
    subprocess.call(cmd)
    
    # now do the SON files
    # get the dynamic files - do the DJFMA files first
    MSLP_SON = output_path + "/ga.pd/field8/" + umid + "_field8_SON.nc"
    WIND_SON = output_path + "/ga.pd/field50/" + umid + "_field50_SON.nc"
    PRECIP_SON = output_path + "/ga.pd/field90/" + umid + "_field90_SON.nc"


    cmd = [EXE_PATH, "-m", MSLP_SON, "-M", MSLP_VNAME,
                     "-w", WIND_SON, "-W", WIND_VNAME,
                     "-p", PRECIP_SON, "-P", PRECIP_VNAME,
                     "-q", POP_FILE, "-Q", POP_VNAME,
                     "-e", REMAP_FILE, "-l", LSM_FILE, "-L", LSM_VNAME,
                     "-t", "7", "-r", "1000.0", "-i",
                     SON_track_name, "-o", EVENT_umid]
    subprocess.call(cmd)

###############################################################################

if __name__ == "__main__":
    opts, args = getopt.getopt(sys.argv[1:], 'y:', ['year='])

    for opt, val in opts:
        if opt in ['--year', '-y']:
            year = val

    # global variables:
    LEV = 7
    EX_LEV = 4
    LS_LEV = 3

    # this is a list of the fields we want to extract:
    # ga.pd indicates daily output from the regional model
    # the fields are:
    # field90 - daily precip - this will be replicated to 6 hourly (4 times replication)
    # field1  - geopotential height (at 700hPa) - 6 hourly instantaneous
    # field50 - 10m wind speed - 6 hour maximum
    # field8  - MSLP - 6 hourly instantaneous
    
    field_list = [["ga.pd", "field90"],
                  ["ga.pd", "field1"],
                  ["ga.pd", "field50"],
                  ["ga.pd", "field8"]]
    
    # create a temporary directory - do we have permission?
    temp_dir = tempfile.mkdtemp()
    
    # set the output directory
    output_dir = "/Users/neil/wah_data/"
    n_valid = 8

    # get the batch list
    batch_file_name = "cpdn_batch/Batch48_1.txt"
    batch_list = read_batch_file(batch_file_name, year)
    for b in batch_list:
        # build the urls for the batch - we want files 1->5 (Dec->Apr)
        # and files 10->12 (Sep->Oct)
        urls = get_urls_for_batch(b)
        for u in urls:
            # extract the data
            #extract(u, field_list, output_dir, temp_dir, n_valid)
            pass
        # concatenate the Dec->Apr files and the Sep->Nov files
        u = urls[0]
        DJFMA_fname, SON_fname, output_path = concatenate(u, field_list, output_dir)
        # do the regridding
        DJFMA_regrid_name, SON_regrid_name = regrid(DJFMA_fname, SON_fname, LEV)
        # find the extrema
        DJFMA_extrema_name, SON_extrema_name = extrema(DJFMA_regrid_name, SON_regrid_name, LEV, EX_LEV, LS_LEV)
        # find the tracks
        DJFMA_track_name, SON_track_name = track(DJFMA_extrema_name, SON_extrema_name)
        # find the events, need the year and umid to form the filename
        umid, boinc, year = get_umid_boinc_year(u)
        events(DJFMA_track_name, SON_track_name, LEV, EX_LEV, LS_LEV, umid, boinc, year)
    
    # clean up the temporary directory
    shutil.rmtree(temp_dir)

