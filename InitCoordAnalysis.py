#!/usr/bin/python

################################################################################
#
# Where it all begins. This script calls the entire data analysis pipeline.
#
# By Chris Ing, 2013 for Python 2.7
#
################################################################################
from CoordAnalysis import *
from argparse import ArgumentParser

def main(args):
    pre_prefix="1st/"
    print "Sorting the ion coordination data"
    data_f = process_input(filenames=args.filenames,
                                          num_cols=args.num_cols,
                                          max_ions=args.max_ions,
                                          remove_frames=args.remove_frames,
                                          traj_col=args.traj_col,
                                          sort_col=args.sort_col,
                                          add_time=args.add_time)

    #write_columns(data_f, outfile="test")

    # Ion occupancy in the channel (governed by the MDAnalysis selection radii)
    print "Channel Occupancy"
    print occ_counter(data_f, num_cols=args.num_cols,
                      traj_col=args.traj_col, sf_col=[],
                      prefix=pre_prefix+"chanocc")

    write_occ_vs_time(data_f, num_cols=args.num_cols,
                      prefix=pre_prefix+"chanocc")

    # Same thing but now we pass the SF column list.
    print "SF Occupancy using columns", args.sf_col
    print occ_counter(data_f, num_cols=args.num_cols,
                      traj_col=args.traj_col, sf_col=args.sf_col,
                      prefix=pre_prefix+"sfocc")

    write_occ_vs_time(data_f, num_cols=args.num_cols,
                      sf_col=args.sf_col, prefix=pre_prefix+"sfocc")

    print "Writing multiple 1D histograms for ion coordination integers"
    write_coord_histograms(data_f, coord_cols=args.coord_cols,
                          num_cols=args.num_cols,
                          sort_col=args.sort_col, prefix=pre_prefix+"1dhisto")

    print "Writing multiple 1D histograms for zero or non-zero coordination"
    write_group_coord_histograms(data_f, sf_col=args.sf_col,
                          num_cols=args.num_cols,
                          sort_col=args.sort_col, prefix=pre_prefix+"1dhisto")

    print "Computing the state stream using regular expressions"
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

    print "Regex Macrostate Occupancy"
    print regex_counter(data_f_padded, data_f_regex,
                        num_ions_map=args.num_ions_map,
                        traj_col=args.traj_col)

    print "Writing 1D histograms for regex coordination labels"
    write_regex_histograms(data_f_padded, data_f_regex,
                           traj_col=args.traj_col,
                           num_cols=args.num_cols,
                           max_ions=args.max_ions,
                           sort_col=args.sort_col,
                           prefix=pre_prefix+"1dhisto")


if __name__ == '__main__':
    parser = ArgumentParser(
    description='This script parses input columnular ASCII data\
    and makes it nice and pretty for subsequent analysis.')

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
    help='a zero inclusive column number to sort your row on, typically x,y,z')
    parser.add_argument(
    '-t', dest='traj_col', type=int, default=11,
    help='a zero inclusive column number that contains the run number')
    parser.add_argument(
    '-o', dest='outfile', type=str, default=None,
    help='the file to output the sorted padding output of all input files')
    parser.add_argument(
    '--addtime', dest='add_time', action="store_true", default=False,
    help='an optional argument to add time columns to each ion grouping')
    parser.add_argument(
    '-coord', dest='coord_cols', type=int, nargs="+", default=[4,5,6,7,8,9,10],
    help='all the columns with coordination counts')
    parser.add_argument(
    '-n', dest='num_ions_map', type=int, nargs="+", default=[1,2,2,3,0],
    help='list of integer ion counts in SF for each regex value + 1 for extra')

    # The following arguments are used for regex state stream processing
    parser.add_argument(
    '-i', dest='regex', type=str, nargs="+",
    help='a list of regex values in quotes for state stream processing')
    parser.add_argument(
    '-sc', dest='sort_cut', type=float, default=0.0,
    help='a value on the sort_col range to classify zero coordinated data')
    parser.add_argument(
    '-sf', dest='sf_col', type=int, nargs="+", default=[5,6],
    help='the coordination integer columns that define the selectivity filter')
    args = parser.parse_args()

    main(args)
