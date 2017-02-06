from datetime import *

#*****************************************************************************

def get_date(units_since, date_ref, calendar_type):

    # get the date from the date_ref
    drs = date_ref.split()
    time_scale = 0.0
    # get a time scalar to convert to hours
    if drs[0] == "days":
        time_scale = 1.0
    elif drs[0] == "hours":
        time_scale = 1.0/24.0
    elif drs[0] == "minutes":
        time_scale = 1.0/(24.0*60.0)

    # get the ref_time
    ref_time_string = drs[2] + ":" + drs[3]
    ref_time = datetime.strptime(ref_time_string, "%Y-%m-%d:%H:%M:%S")

    if calendar_type in ["standard", "gregorian"]:
        tgt_time = ref_time + timedelta(days = units_since * time_scale)
        return (tgt_time.year, tgt_time.month, tgt_time.day, tgt_time.hour)
    elif calendar_type in ["360_day", "360day"]:
    	ref_days = (ref_time.year) * 360 + (ref_time.month-1) * 30 + (ref_time.day-1) + (ref_time.hour) * 1.0/24
        n_evt_days = units_since * time_scale + ref_days
        evt_yr = int(n_evt_days / 360.0)
        n_evt_days -= evt_yr * 360
        evt_mon = int(n_evt_days / 30.0)
        n_evt_days -= evt_mon * 30
        evt_day = int(n_evt_days)
        n_evt_days -= evt_day
        evt_hour = int(n_evt_days * 24.0)
        return (evt_yr, evt_mon+1, evt_day+1, evt_hour)
        

#****************************************************************************

if __name__ == "__main__":
	# some tests for values I already know
	print get_date(119.5417, "days since 1989-12-1 00:00:00", "standard")
	print get_date(9360.0, "days since 1959-12-1 00:00:00", "360_day")
	print get_date(9630.25, "days since 1959-12-1 00:00:00", "360_day")