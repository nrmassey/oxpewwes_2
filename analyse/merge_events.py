#! /usr/bin/env python

###############################################################################
# Program: merge_events.py
# Purpose: Merge two events in an event set - useful to match historical storms
#          to those in XWS where the automated storm tracking algorithm may have
#          failed & created two storm tracks, rather than just one
# Author : Neil Massey
# Date   : 04/04/2016
# Contact: neil.massey@ouce.ox.ac.uk
###############################################################################

import sys, getopt
from read_write_event import read_event, write_event
from haversine import haversine

def list_overlapping_events(input_fname, sr):
    """List the overlapping events in the event set.
       Information will include the indices of the two overlapping events and the timesteps at which they overlap.
       This information can then be used to perform the merge.
    """
    # load the event
    event = read_event(input_fname)

    # loop over all the events / tracks
    for tr_A in range(0, event.n_evts):
        # get the start and end timestep of track A
        tr_A_st = event.track_time[tr_A][0]
        tr_A_ed = event.track_time[tr_A][-1]
        # loop over all the other tracks
        for tr_B in range(tr_A+1, event.n_evts):
            tr_B_st = event.track_time[tr_B][0]
            tr_B_ed = event.track_time[tr_B][-1]

            # This bit of code mirrors that in tracker.cpp in the main tracker app
            # four overlapping scenarios handled by three clauses:
            # +----+     +----+      +----+       +-+          tr_A
            #   +----+     +----+     +--+      +-----+        tr_B
            if ((tr_A_ed >= tr_B_st and tr_A_ed <= tr_B_ed) or
                (tr_A_st >= tr_B_st and tr_A_st <= tr_B_ed) or
                (tr_A_st >= tr_B_st and tr_A_ed <= tr_B_ed)):

                # tracks must have one point at one timestep that is 
                # < search radius between each track
                within_radius = False
                for b in range(0, len(event.track_time[tr_B])):
                    # get distance between this point and the last track point of A
                    # and then the first track point of A
                    dist_AfB = haversine(event.track_longitude[tr_A][0], event.track_latitude[tr_A][0], 
                                         event.track_longitude[tr_B][b], event.track_latitude[tr_B][b])
                    dist_AlB = haversine(event.track_longitude[tr_A][-1], event.track_latitude[tr_A][-1],
                                         event.track_longitude[tr_B][b], event.track_latitude[tr_B][-1])
                    if (dist_AfB < sr or dist_AlB < sr):
                        within_radius = True
                        break
                if (within_radius):
                    tr_A_times = ",".join([str(event.track_time[tr_A][x]) for x in range(0, event.track_time[tr_A].shape[0])])
                    tr_B_times = ",".join([str(event.track_time[tr_B][x]) for x in range(0, event.track_time[tr_B].shape[0])])
                    print "Overlapping tracks: " + str(tr_A) + "[" + tr_A_times + "] " + str(tr_B) + "[" + tr_B_times + "]"


if __name__ == "__main__":
    """Merge two events by merging the tracks and associated fields.
       The merged track will consist of subtracks from two existing tracks: m and n.
       For track m, the timesteps included will be 0->t
       For track n, the timesteps included will be s->len(n)
       Use -l to list which tracks can be merged and at which timesteps
        Options are:
        -i: input file
        -m: index of first event to merge
        -t: timestep of first event to truncate to
        -n: index of second event to merge
        -s: timestep of second event to begin track from
        -l: list possible merged tracks
    """
    opts, args = getopt.getopt(sys.argv[1:], 'i:m:n:s:t:l', ['input=', 'event1=', 'event2=', 'tstep1=', 'tstep2='])

    input_file = ""
    search_radius = 100.0
    first_event_idx = -1
    second_event_idx = -1
    first_event_tstep = -1
    second_event_tstep = -1
    list_events = True
 
    for opt, val in opts:
        if opt in ['--input', '-i']:
            input_file = val
        if opt in ['--event1', '-m']:
            first_event_idx = int(val)
        if opt in ['--event2', '-n']:
            second_event_idx = int(val)
        if opt in ['--tstep1', '-t']:
            first_event_tstep = int(val)
        if opt in ['--tstep2', '-s']:
            second_event_tstep = int(val)
        if opt in ['--sr', '-r']:
            search_radius = float(val)
        if opt in ['--list', '-l']:
            list_events = True

    if input_file == "":
        print "Missing input file (input|i)"
        sys.exit()

    if list_events:
        list_overlapping_events(input_file, search_radius)
    else:
        exit = False
        if first_event_idx == -1:
            print "Missing first event option (event1|m)"
            exit = True
        if first_event_idx == -1:
            print "Missing second event option (event2|n)"
            exit = True
        if first_event_idx == -1:
            print "Missing first timestep option (tstep1|t)"
            exit = True
        if first_event_idx == -1:
            print "Missing second timestep option (tstep2|s)"
            exit = True
        if exit:
            sys.exit()
 
