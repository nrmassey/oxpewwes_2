import numpy

###############################################################################

def load_data(nc_fh, var_name, index):
    # load data that has been packed to a byte or short and scaled
    # NB in netCDF4 python library the scaling and offset are performed automagically
    var = nc_fh.variables[var_name]
    # get the missing value
    mv = var._FillValue
    # get the data and mask
    D = numpy.ma.masked_equal(var[index], mv)
    return D
