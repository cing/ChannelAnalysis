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

    #write_columns(data_f, outfile=pre_prefix+"test")

    # Ion occupancy in the channel (governed by the MDAnalysis selection radii)
    print "Channel Occupancy"
    print occ_counter(data_f, num_cols=args.num_cols,
                      traj_col=args.traj_col, sf_col=[],
                      prefix=pre_prefix+"chanocc")

    write_occ_vs_time(data_f, num_cols=args.num_cols,
                      traj_col=args.traj_col, prefix=pre_prefix+"chanocc")

    # Same thing but now we pass the SF column list.
    print "SF Occupancy using columns", args.sf_col
    print occ_counter(data_f, num_cols=args.num_cols,
                      traj_col=args.traj_col, sf_col=args.sf_col,
                      prefix=pre_prefix+"sfocc")

    write_occ_vs_time(data_f, num_cols=args.num_cols, traj_col=args.traj_col,
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

    print "State Transition Counting"
    data_state_trans = state_transitions(data_f_padded, data_f_regex,
                            time_col=args.time_col,
                            time_increment=args.time_increment,
                            traj_col=args.traj_col)

    print "Computing Regex State Occupancies for macrostate graph"
    data_f_occupancy = regex_counter(data_f_padded, data_f_regex,
                        num_ions_map=args.num_ions_map,
                        traj_col=args.traj_col)
    # We want to extract only the mean occupancies in percent (2nd last entry)
    data_f_mean_pop = data_f_occupancy[-2][1:]

    print "State Transition Graph Building and Writing"
    #print data_state_trans
    state_draw_map = [(0.0, 0),(1.0, 0),(2.0, 0),(3.0, 0),
                          (0.5,-1),(1.5,-1),(2.5,-1),
                               (1.0,-2),(2.0,-2),(3.0,-2),
                                   (1.5,-3),(2.5,-3),
                                        (2.0,-4)]

    data_state_trans_graph = build_macrostate_graph(data_state_trans,
                                                    state_draw_map,
                                                    pop_map=data_f_mean_pop)

    print "Macrostate Graph Writing"
    write_macrostate_graph(data_state_trans_graph,
                            outfile=pre_prefix+"graph_regex_macro.pdf")

    print "Microstate Graph Writing"
    write_microstate_graph(data_f_padded, data_f_regex,
                           time_col=args.time_col,
                           time_increment=args.time_increment,
                           traj_col=args.traj_col,
                           outfile=pre_prefix+"graph_regex_micro.gexf")

    print "State Transition w/ Intermediate Counting"
    print state_intermediate_transitions(data_f_padded, data_f_regex,
                                         time_col=args.time_col,
                                         time_increment=args.time_increment,
                                         traj_col=args.traj_col)


if __name__ == '__main__':
    parser = ArgumentParser(
    description='This script parses input columnular ASCII data\
    and makes it nice and pretty for subsequent analysis.')

    # These arguments are required for input processing and are required
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
    '--addtime', dest='add_time', action="store_true", default=False,
    help='an optional argument to add time columns to each ion grouping')

    # The following is used for 1D coordination histograms
    parser.add_argument(
    '-coord', dest='coord_cols', type=int, nargs="+", default=[4,5,6,7,8,9,10],
    help='all the columns with coordination counts')

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
    parser.add_argument(
    '-n', dest='num_ions_map', type=int, nargs="+", default=[1,2,2,3,0],
    help='list of integer ion counts in SF for each regex value + 1 for extra')

    # The following arguments are used for regex state transition processing
    parser.add_argument(
    '-timecol', dest='time_col', type=int, default=0,
    help='a zero inclusive column number with the frame/timestep number')
    parser.add_argument(
    '-dt', dest='time_increment', type=int, default=1,
    help='the difference in the time column between steps')

    args = parser.parse_args()

    main(args)
