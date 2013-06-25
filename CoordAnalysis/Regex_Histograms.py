#!/usr/bin/python
###############################################################################
#
# This script extracts coordination information as a function of a cartesian
# coordination column and creates 1D histograms as output. It combines
# the previous script Extract_Lines_with_Coordination_Integer.py and 1DHisto.py
# as well as all the E_only_Coordination_SF.py-type scripts.
#
# Example: For 13/26/39/...-column data with type like:
#          1.0 -0.13 -0.193 0.522 0.0 1.0 0.0 0.0 0.0 2.0 9.0 2.0 1748.0
#
#          The following command would remove 2000 lines from the input
#          and produce a large number of 1D Histograms to be plotted
#          in an external plotting program (split by number of ions):
#          python Regex_Histograms.py -f f1.out -m 3 -c 13 -remove 2000
#                                  -i "(.[^-0+]|[^-0+].)[-0][-0][-0][-0]"
#                                  "\+\+[-0][-0][-0][-0]"
#                                  "(.[^-+0]|[^-+0].)(.[^-+0]|[^-+0].)[-0][-0]"
#                                  "\+\+(.[^-0+]|[^-0+].)[-0][-0]"
#                                  -sf 5 6
#
# By Chris Ing, 2013 for Python 2.7
#
###############################################################################
from argparse import ArgumentParser
from collections import defaultdict
from numpy import histogram
from itertools import product
from re import match
from Ion_Preprocessor import *

def compute_regex_histograms(data_floats, data_regex, num_cols=13, pad_col=4,
                           sort_col=3, max_ions=3, traj_col=11, histmin=-1.50,
                           histmax=1.5, histbins=300, prefix=None):

    # This is an epic datatype with the 1st key as the regex string
    # of interest, the 2nd key is the ion number within that grouping
    # and the value is a list of sort values where that integer was observed.
    coord_sortvals=defaultdict(lambda: defaultdict(list))

    # These are dictionaries of dictionaries where the key is a regex id
    # and the list is a axial probability or associated z value.
    hist_per_regex = defaultdict(list)
    z_per_regex = defaultdict(list)

    # the data_regex datatype is a list of tuples: (state_label, regex_int)
    data_regex_ids = [data[1] for data in data_regex]

    for line, regex_state in zip(data_floats, data_regex_ids):
        traj_id = line[traj_col]

        for ion_index, ion in enumerate(chunker(line,num_cols)):
            if ion[pad_col] != "-":
                coord_sortvals[regex_state][ion_index].append(ion[sort_col])

    # Iterate over the data structure that was built and output
    # all the necessary histograms.
    for coord_col, coord_dict in coord_sortvals.iteritems():

        combined_sort_vals = []
        for coord_int, sort_vals in coord_dict.iteritems():

            combined_sort_vals.extend(sort_vals)
            histo, edges = histogram(sort_vals, range=[histmin, histmax],
                                     bins=histbins, normed=False)

            if prefix != None:
                # Print out histograms separately for each ion
                with open(prefix+"_regex"+str(coord_col)+
                                 "_ion"+str(coord_int),"w") as out:
                    for xval, yval in zip(edges,histo):
                        out.write(str(xval)+" "+str(yval)+"\n")

            hist_per_regex[str(coord_col)+
                               "_"+str(coord_int)].append(histo)
            z_per_regex[str(coord_col)+"_"+str(coord_int)].append(edges)

        # Compute a histogram for all ions to act as an envelope
        histo, edges = histogram(combined_sort_vals, range=[histmin, histmax],
                                     bins=histbins, normed=False)

        if prefix != None:
            with open(prefix+"_regex"+str(coord_col)+
                      "_ionALL","w") as out:
                for xval, yval in zip(edges,histo):
                    out.write(str(xval)+" "+str(yval)+"\n")

        hist_per_regex[str(coord_col)+"_ALL"].append(histo)
        z_per_regex[str(coord_col)+"_ALL"].append(edges)

    return (hist_per_regex.items(), z_per_regex.items())

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
    '-sc', dest='sort_cut', type=float, default=0.0,
    help='a value on the sort_col range to classify zero coordinated data')
    parser.add_argument(
    '-t', dest='traj_col', type=int, default=11,
    help='a zero inclusive column number that contains the run number')
    parser.add_argument(
    '-sf', dest='sf_col', type=int, nargs="+", default=[5,6],
    help='the coordination integer columns that define the selectivity filter')
    parser.add_argument(
    '-p', dest='prefix', type=str, default="1dhisto",
    help='this is the prefix of the histogram output')
    parser.add_argument(
    '--addtime', dest='add_time', action="store_true", default=False,
    help='an optional argument to add time columns to each ion grouping')
    parser.add_argument(
    '-i', dest='regex', type=str, nargs="+", required=True,
    help='a list of regex values in quotes')
    args = parser.parse_args()

    data_f_padded = process_input(filenames=args.filenames,
                                          num_cols=args.num_cols,
                                          max_ions=args.max_ions,
                                          remove_frames=args.remove_frames,
                                          traj_col=args.traj_col,
                                          sort_col=args.sort_col,
                                          add_time=args.add_time,
                                          padded=True)

    data_f_regex = regex_columns(data_f_padded, regex_strings=args.regex,
                                num_cols=args.num_cols,
                                sort_col=args.sort_col,
                                sort_cut=args.sort_cut,
                                sf_col=args.sf_col,
                                max_ions=args.max_ions)

    print "Writing 1D histograms for regex coordination labels"
    print compute_regex_histograms(data_f_padded, data_f_regex,
                                   traj_col=args.traj_col,
                                   num_cols=args.num_cols,
                                   max_ions=args.max_ions,
                                   sort_col=args.sort_col,
                                   prefix=args.prefix)
