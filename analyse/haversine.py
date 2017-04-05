import math

EARTH_R = 6371 # kilometers

def haversine(lon1, lat1, lon2, lat2, R=EARTH_R):
    """Function to calculate the haversine distance between two lat / lon points"""
    # quick check to prevent very small distances at the poles
    if (lat1 == 90.0):
        lat1 = 90.5
    if (lat1 == -90.0):
        lat1 = -90.5
    if (lat2 == 90.0):
        lat2 = 90.5
    if (lat2 == -90.0):
        lat2 = -90.5

    # convert lon / lat into radians
    lon1_r = math.radians(lon1)
    lat1_r = math.radians(lat1)
    lon2_r = math.radians(lon2)
    lat2_r = math.radians(lat2)

    # differences
    d_lon = lon2_r - lon1_r
    d_lat = lat2_r - lat1_r

    # component parts of equation
    a = math.sin(d_lat/2.0)
    b = math.sin(d_lon/2.0)
    c = a*a + math.cos(lat1_r) * math.cos(lat2_r) * b*b
    # don't scupper the arcsin
    if (c > 1.0):
        c = 1.0
    d = R * 2 * math.asin(math.sqrt(c))
    return d
