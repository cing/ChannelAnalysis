#!/usr/bin/python
###############################################################################
#
# This script extracts coordination information as a function of a cartesian
# coordination column and performs a 1D histogram of the output. It combines
# the previous script Extract_Lines_with_Coordination_Integer.py and 1DHisto.py
# as well as all the E_only_Coordination_SF.py-type scripts.
#
# Example: For 13/26/39/...-column data with type like:
#          1.0 -0.13 -0.193 0.522 0.0 1.0 0.0 0.0 0.0 2.0 9.0 2.0 1748.0
#
#          The following command would remove 2000 lines from the input
#          and produce a large number of 1D Histograms to be plotted
#          in an external plotting program:
#          python Coordination_Histograms.py -f f1.out -m 3 -c 13 -remove 2000
#                                         -coord 4 5 6 7 8 9 10
#
# By Chris Ing, 2013 for Python 2.7
#
###############################################################################
from argparse import ArgumentParser
from collections import defaultdict
from numpy import histogram, histogram2d, sqrt, linspace, zeros, digitize, log
from numpy import where, isinf
from itertools import combinations
from Ion_Preprocessor import *

# This returns the sort_column as a time series, useful
# for making scatterplot time series of ionic positions.
def compute_ion_timeseries(data_lines, sort_col, traj_col,
                           pad_col=4, num_cols=13):

    # These are dictionaries of dict where the key is the traj_number
    # and the subdict is ion_number and te value is a LIST of ion positions,
    # or associated time values in the case of the associated time_per_traj
    ion_pos_per_traj = defaultdict(dict)
    time_per_traj = defaultdict(dict)

    for line in data_lines:
        for ion_num, ion in enumerate(list(chunker(line,num_cols))):
            # In the case that the ion grouping has a hypen
            # we know that's a padded column and must be excluded.
            if ion[pad_col] != "-":
                sort_val = ion[sort_col]
                traj_id = ion[traj_col]
                if ion_num not in ion_pos_per_traj[traj_id]:
                    ion_pos_per_traj[traj_id][ion_num] = [sort_val]
                    time_per_traj[traj_id][ion_num] = [line[0]]
                else:
                    ion_pos_per_traj[traj_id][ion_num].append(sort_val)
                    time_per_traj[traj_id][ion_num].append(line[0])

    return (dict(ion_pos_per_traj), dict(time_per_traj))

# To accompany compute_ion_timeseries, this function computes ion
# position histograms per ion without taking into consideration
# the coordination or macrostate that ion is occupying.
def compute_position_histograms(data_lines, sort_col, traj_col,
                                pad_col=4, num_cols=13,
                                histmin=-1.50, histmax=1.5,
                                histbins=300, histcount_cut=200):

    # Power datatypes son. The structure is: traj_id -> ion_num -> z_vals
    ion_sortvals = defaultdict(lambda: defaultdict(list))

    # These are dictionaries of lists where the key is a coord_col
    # and the list is a axial probability or associated z value.
    coord_hist_per_ion = defaultdict(list)
    z_per_ion= defaultdict(list)

    for line in data_lines:
        for ion_num, ion in enumerate(list(chunker(line,num_cols))):
            # In the case that the ion grouping has a hypen
            # we know that's a padded column and must be excluded.
            if ion[pad_col] != "-":
                sort_val = ion[sort_col]
                traj_id = ion[traj_col]
                ion_sortvals[traj_id][ion_num].append(sort_val)

    for traj_id, ion_dict in ion_sortvals.iteritems():
        for ion_num, ion_vals in ion_dict.iteritems():
            if len(ion_vals) > histcount_cut:
                histo, edges = histogram(ion_vals, range=[histmin, histmax],
                                     bins=histbins, normed=True)
                coord_hist_per_ion[traj_id].append(histo)
                z_per_ion[traj_id].append(edges)

    return (dict(coord_hist_per_ion), dict(z_per_ion))

# Exactly as above, except we don't histogram each ion individually,
# hence the confusing function name change of "histograms" to "histogram".
def compute_position_histogram(data_lines, sort_col, traj_col,
                                pad_col=4, num_cols=13,
                                histmin=-1.50, histmax=1.5,
                                histbins=300,):

    # Power datatypes son. The structure is: traj_id -> ion_num -> z_vals
    ion_sortvals = defaultdict(list)

    # These are dictionaries of lists where the key is a coord_col
    # and the list is a axial probability or associated z value.
    coord_hist_per_ion = defaultdict(list)
    z_per_ion= defaultdict(list)

    for line in data_lines:
        for ion in chunker(line,num_cols):
            # In the case that the ion grouping has a hypen
            # we know that's a padded column and must be excluded.
            if ion[pad_col] != "-":
                sort_val = ion[sort_col]
                traj_id = ion[traj_col]
                ion_sortvals[traj_id].append(sort_val)

    for traj_id, ion_vals in ion_sortvals.iteritems():
        histo, edges = histogram(ion_vals, range=[histmin, histmax],
                                bins=histbins, normed=True)
        coord_hist_per_ion[traj_id].extend(list(histo))
        z_per_ion[traj_id].extend(list(edges))

    return (dict(coord_hist_per_ion), dict(z_per_ion))

# Yet another way of making a histogram, we don't look at the traj_col
# at all, we simply add all the ionic data together and do 1 big boy.
def compute_allion_histogram(data_lines, sort_col,
                                pad_col=4, num_cols=13,
                                histmin=-1.50, histmax=1.5,
                                histbins=300, normed=True, prefix=None):

    # Power datatypes son. The structure is: traj_id -> ion_num -> z_vals
    ion_sortvals = []

    # These are dictionaries of lists where the key is a coord_col
    # and the list is a axial probability or associated z value.
    coord_hist_per_ion = defaultdict(list)
    z_per_ion= defaultdict(list)

    for line in data_lines:
        for ion in chunker(line,num_cols):
            # In the case that the ion grouping has a hypen
            # we know that's a padded column and must be excluded.
            if ion[pad_col] != "-":
                sort_val = ion[sort_col]
                ion_sortvals.append(sort_val)

    histo, edges = histogram(ion_sortvals, range=[histmin, histmax],
                            bins=histbins, normed=normed)

    if prefix != None:
        with open(prefix+"_allion","w") as out:
            for xval, yval in zip(edges,histo):
                out.write(str(xval)+" "+str(yval)+"\n")

    coord_hist_per_ion["ALL"].extend(list(histo))
    z_per_ion["ALL"].extend(list(edges))

    return (dict(coord_hist_per_ion), dict(z_per_ion))

# This function exists to compute the square root of the sum of two squared
# values in each line. For practical purposes, this is useful for
# make a time series of r (sqrt(x^2+y^2)) and returning it in the same
# format as the normal compute_ion_timeseries
def compute_ion_sqrt_timeseries(data_lines, square_cols, traj_col,
                           pad_col=4, num_cols=13):

    # These are dictionaries of dict where the key is the traj_number
    # and the subdict is ion_number and te value is a LIST of ion positions,
    # or associated time values in the case of the associated time_per_traj
    ion_pos_per_traj = defaultdict(dict)
    time_per_traj = defaultdict(dict)

    for line in data_lines:
        for ion_num, ion in enumerate(list(chunker(line,num_cols))):
            # In the case that the ion grouping has a hypen
            # we know that's a padded column and must be excluded.
            if ion[pad_col] != "-":
                square_vals = [ion[sq_col]**2 for sq_col in square_cols]
                root_val = sqrt(sum(square_vals))
                traj_id = ion[traj_col]
                if ion_num not in ion_pos_per_traj[traj_id]:
                    ion_pos_per_traj[traj_id][ion_num] = [root_val]
                    time_per_traj[traj_id][ion_num] = [line[0]]
                else:
                    ion_pos_per_traj[traj_id][ion_num].append(root_val)
                    time_per_traj[traj_id][ion_num].append(line[0])

    return (dict(ion_pos_per_traj), dict(time_per_traj))

# To accompany compute_ion_sqrt_timeseries, this function computes ion
# position histograms per ion without taking into consideration
# the coordination or macrostate that ion is occupying.
def compute_position_sqrt_histograms(data_lines, square_cols, traj_col,
                                pad_col=4, num_cols=13,
                                histmin=0, histmax=1.0,
                                histbins=300, histcount_cut=200):

    # Power datatypes son. The structure is: traj_id -> ion_num -> z_vals
    ion_sortvals = defaultdict(lambda: defaultdict(list))

    # These are dictionaries of lists where the key is a coord_col
    # and the list is a axial probability or associated z value.
    coord_hist_per_ion = defaultdict(list)
    z_per_ion= defaultdict(list)

    for line in data_lines:
        for ion_num, ion in enumerate(list(chunker(line,num_cols))):
            # In the case that the ion grouping has a hypen
            # we know that's a padded column and must be excluded.
            if ion[pad_col] != "-":
                square_vals = [ion[sq_col]**2 for sq_col in square_cols]
                root_val = sqrt(sum(square_vals))
                traj_id = ion[traj_col]
                ion_sortvals[traj_id][ion_num].append(root_val)

    for traj_id, ion_dict in ion_sortvals.iteritems():
        for ion_num, ion_vals in ion_dict.iteritems():
            if len(ion_vals) > histcount_cut:
                histo, edges = histogram(ion_vals, range=[histmin, histmax],
                                     bins=histbins, normed=True)
                coord_hist_per_ion[traj_id].append(histo)
                z_per_ion[traj_id].append(edges)

    return (dict(coord_hist_per_ion), dict(z_per_ion))

# this function looks at all ions at a timestep and checks what
# coordination they have, grouping them into a list for generating
# multiple 1D Histograms. coord_cols is a list of columns that contain
# coordination counts.
def compute_coord_histograms(data_lines, coord_cols=[4,5,6,7,8,9,10],
                             num_cols=13,
                             sort_col=3, pad_col=4, histmin=-1.50, histmax=1.5,
                             histbins=300, prefix=None):

    # This is an epic datatype with the 1st key as the coordination column
    # of interest, the 2nd key is the integer coordination of interest
    # and the value is a list of values where that integer was observed.
    coord_sortvals = defaultdict(lambda: defaultdict(list))

    # These are dictionaries of lists where the key is a coord_col
    # and the list is a axial probability or associated z value.
    coord_hist_per_col = defaultdict(list)
    z_per_col = defaultdict(list)

    for line in data_lines:
        for ion in chunker(line,num_cols):
            # In the case that the ion grouping has a hypen
            # we know that's a padded column and must be excluded.
            if ion[pad_col] != "-":
                sort_val = ion[sort_col]

                # He we accumulate all the sort_vals based on integer
                # coordination into list format. This will be used
                # to perform 1D Histogramming.
                for coord_col in coord_cols:
                    coord_sortvals[coord_col][ion[coord_col]].append(sort_val)

    # Iterate over the data structure that was built and output
    # all the necessary histograms.
    for coord_col, coord_dict in coord_sortvals.iteritems():
        for coord_int, sort_vals in coord_dict.iteritems():
            histo, edges = histogram(sort_vals, range=[histmin, histmax],
                                     bins=histbins, normed=False)

            if prefix != None:
                with open(prefix+"_coordcol"+str(coord_col)+
                                 "_coordint"+str(coord_int),"w") as out:
                    for xval, yval in zip(edges,histo):
                        out.write(str(xval)+" "+str(yval)+"\n")

            coord_hist_per_col[str(coord_col)+
                               "_"+str(coord_int)].append(histo)
            z_per_col[str(coord_col)+"_"+str(coord_int)].append(edges)

    return (coord_hist_per_col.items(), z_per_col.items())

# similar to above, but this script groups the data based on zero
# and non-zero values of the passed coord_cols. For example,
# by passing coord_cols=[5,6], this script will output the following files:
# 1) 5: >0, 6: 0    (col5 coordination only) as file suffix 10
# 2) 5: >0, 6: >0   (col5 and col6 coordination) as file suffix 11
# 3) 5: 0,  6: >0   (col6 coordination only) as file suffix 01
# 4) 5: 0,  6: 0    (col5 and col6 no coordination) as file suffice 00
# 5) -----------    (all ions with coordination, sum of 1,2,3 above) as ++
def compute_group_coord_histograms(data_lines, sf_col=[5,6],
                                   num_cols=13,
                                   sort_col=3, pad_col=4, histmin=-1.50,
                                   histmax=1.5, histbins=300, prefix=None):

    # This is a datatype where the 1st key is the coordination group id
    # and the value is a list of values where that group id was observed.
    coord_sortvals = defaultdict(list)

    # These are dictionaries of dictionaries where the key is a coord_col
    # and the list is a axial probability or associated z value.
    coord_hist_per_col = defaultdict(list)
    z_per_col = defaultdict(list)

    for line in data_lines:
        for ion in chunker(line,num_cols):
            # In the case that the ion grouping has a hypen
            # we know that's a padded column and must be excluded.
            if ion[pad_col] != "-":
                sort_val = ion[sort_col]

                # Extract the coordination using the columns passed
                coords = [ion[col] for col in sf_col]
                # The same as above but in string form with booleans
                coords_str = "".join([str(int(coord>0)) for coord in coords])

                coord_sortvals[coords_str].append(sort_val)
                if any([coord>0 for coord in coords]):
                    coord_sortvals['+'*len(coords)].append(sort_val)

    # Iterate over the data structure that was built and output
    # all the necessary histograms.
    for group_id, sort_vals in coord_sortvals.iteritems():
        histo, edges = histogram(sort_vals, range=[histmin, histmax],
                                 bins=histbins, normed=False)

        if prefix != None:
            with open(prefix+"_groupcoord"+str(group_id),"w") as out:
                for xval, yval in zip(edges,histo):
                    out.write(str(xval)+" "+str(yval)+"\n")

        coord_hist_per_col[str(group_id)].append(histo)
        z_per_col[str(group_id)].append(edges)

    return (coord_hist_per_col.items(), z_per_col.items())

# Returns a histogram for each ion in each of the channel
# ion occupancy states across the whole dataset.
def compute_ionsplit_histograms(data_lines, sort_col,
                                pad_col=4, num_cols=13,
                                histmin=-1.50, histmax=1.5,
                                histbins=300, histcount_cut=200,
                                prefix=None):

    # Power datatypes son. The structure is: ion_occ -> ion_num -> z_vals
    ion_sortvals = defaultdict(lambda: defaultdict(list))

    # These are dictionaries of lists where the key is a coord_col
    # and the list is a axial probability or associated z value.
    coord_hist_per_ion = defaultdict(list)
    z_per_ion= defaultdict(list)

    for line in data_lines:

        temp_ion_count = 0
        for ion in chunker(line,num_cols):
            # In the case that the ion grouping has a hypen
            # we know that's a padded column and must be excluded.
            if ion[pad_col] != "-":
                temp_ion_count += 1

        for ion_num, ion in enumerate(list(chunker(line,num_cols))):
            # In the case that the ion grouping has a hypen
            # we know that's a padded column and must be excluded.
            if ion[pad_col] != "-":
                sort_val = ion[sort_col]
                ion_sortvals[temp_ion_count][ion_num].append(sort_val)

    for occ_id, ion_dict in ion_sortvals.iteritems():
        for ion_num, ion_vals in ion_dict.iteritems():
            if len(ion_vals) > histcount_cut:
                histo, edges = histogram(ion_vals, range=[histmin, histmax],
                                     bins=histbins, normed=True)

                if prefix != None:
                    with open(prefix+"_ionsplit"+str(occ_id)
                              + str(ion_num),"w") as out:

                        for xval, yval in zip(edges,histo):
                            out.write(str(xval)+" "+str(yval)+"\n")

                coord_hist_per_ion[occ_id].append(histo)
                z_per_ion[occ_id].append(edges)

    return (dict(coord_hist_per_ion), dict(z_per_ion))

# Returns a histogram for each ion in each of the channel
# ion occupancy states across the whole dataset.
def compute_ionsplit_2dhistograms(data_lines, sort_col,
                                  pad_col=4, num_cols=13,
                                  occ_cutoff=3, kBT=0.596,
                                  histmin=-1.50, histmax=1.5,
                                  histbins=300,
                                  prefix=None):

    # Power datatypes son. The structure is: ion_occ -> ion_num -> z_vals
    ion_sortvals = defaultdict(lambda: defaultdict(list))

    # These are dictionaries of lists where the key is a coord_col
    # and the list is a axial probability or associated z value.
    coord_2dhist_per_pair = defaultdict(list)
    xedges_per_pair = defaultdict(list)
    yedges_per_pair = defaultdict(list)

    for line in data_lines:

        temp_ion_count = 0
        for ion in chunker(line,num_cols):
            # In the case that the ion grouping has a hypen
            # we know that's a padded column and must be excluded.
            if ion[pad_col] != "-":
                temp_ion_count += 1

        # 2d histograms do not exist for 0 and 1 ion occupancy values.
        if temp_ion_count > 1:

            # Ion occupancy values larger than occ_cutoff are lumped together
            temp_ion_count = min(temp_ion_count, occ_cutoff)

            for ion_num, ion in enumerate(list(chunker(line,num_cols))):
                # In the case that the ion grouping has a hypen
                # we know that's a padded column and must be excluded.
                # Keep in mind that if we had large ion occupancy
                # we won't be plotting the 2d histograms passed the occ_cutoff
                # so it's pointless to generate that data, hence the
                # occ_cutoff check below.
                if ion[pad_col] != "-" and ion_num < occ_cutoff:
                    sort_val = ion[sort_col]
                    ion_sortvals[temp_ion_count][ion_num].append(sort_val)

    for occ_id, ion_dict in ion_sortvals.iteritems():

        for occ_pair in combinations(range(occ_id),2):
            histo, xedges, yedges = histogram2d(ion_dict[occ_pair[0]],
                                                ion_dict[occ_pair[1]],
                                                range=[[histmin, histmax],
                                                       [histmin, histmax]],
                                                bins=[histbins,histbins],
                                                normed=True)

            # Since we're taking the log, remove all 0.0 values.
            low_val_indices = histo <= 0.0
            high_val_indices = histo > 0.0
            histo[low_val_indices] = 5.0
            histo[high_val_indices] = -kBT*log(histo[high_val_indices])

            # Everything must be shifted such that 0 is the true minimum value.
            min_val = min(histo[high_val_indices])
            histo[high_val_indices] = histo[high_val_indices] + abs(min_val)

            if prefix != None:
                with open(prefix + "_ionsplit_occ" + str(occ_id),"w") as out:
                    for xval, yval, zval in zip(xedges, yedges, histo):
                        out.write(str(xval) + " " +
                                  str(yval) + " " + str(zval) + "\n")

            #import pdb
            #pdb.set_trace()
            coord_2dhist_per_pair[occ_id].append(histo)
            xedges_per_pair[occ_id].append(xedges)
            yedges_per_pair[occ_id].append(yedges)

    return (dict(coord_2dhist_per_pair),
            dict(xedges_per_pair),
            dict(yedges_per_pair))

if __name__ == '__main__':
    parser = ArgumentParser(
    description='This script extracts ions with specific integer coordination\
    and computes 1D Histograms along a cartesian axis for those ions')

    parser.add_argument(
    '-f', dest='filenames', type=str, nargs="+", required=True,
    help='a filename of coordination data from MDAnalysis trajectory data')
    parser.add_argument(
    '-m', dest='max_ions', type=int, required=True,
    help='the maximum number of ions in the channel to consider')
    parser.add_argument(
    '-c', dest='num_cols', type=int, default=13,
    help='the number of columns per ion in the input')
    parser.add_argument(
    '-remove', dest='remove_frames', type=int, default=0,
    help='this is a number of frames to remove from the start of the data')
    parser.add_argument(
    '-s', dest='sort_col', type=int, default=3,
    help='a zero inclusive column number to sort your row on')
    parser.add_argument(
    '-t', dest='traj_col', type=int, default=11,
    help='a zero inclusive column number that contains the run number')
    parser.add_argument(
    '-coord', dest='coord_cols', type=int, nargs="+", default=[4,5,6,7,8,9,10],
    help='all the columns with coordination counts')
    parser.add_argument(
    '-sf', dest='sf_col', type=int, nargs="+", default=[5,6],
    help='the coordination integer columns that define the selectivity filter')
    parser.add_argument(
    '-p', dest='prefix', type=str, default="1dhisto",
    help='this is the prefix of the histogram output')
    parser.add_argument(
    '--addtime', dest='add_time', action="store_true", default=False,
    help='an optional argument to add time columns to each ion grouping')
    args = parser.parse_args()

    data_f = process_input(filenames=args.filenames,
                                          num_cols=args.num_cols,
                                          max_ions=args.max_ions,
                                          remove_frames=args.remove_frames,
                                          traj_col=args.traj_col,
                                          sort_col=args.sort_col,
                                          add_time=args.add_time)

    '''
    print "Writing multiple 1D histograms for ion coordination integers"
    print compute_coord_histograms(data_f, coord_cols=args.coord_cols,
                             num_cols=args.num_cols,
                             sort_col=args.sort_col, prefix=args.prefix)

    print "Writing multiple 1D histograms for zero or non-zero coordination"
    print compute_group_coord_histograms(data_f, sf_col=args.sf_col,
                                   num_cols=args.num_cols,
                                   sort_col=args.sort_col, prefix=args.prefix)
    '''

    print "Writing average coordination per z bin arrays"
    compute_avg_coord_histograms(data_f, coord_cols=args.coord_cols,
                                     num_cols=args.num_cols,
                                     sort_col=args.sort_col,
                                     prefix=args.prefix)

