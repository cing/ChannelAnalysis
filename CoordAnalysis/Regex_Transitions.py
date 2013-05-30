#!/usr/bin/python

################################################################################
#
# This script parses a state stream and detects transitions that
# satisfy a particular regular expression for the state in start,
# intermediate, and end states. It allows short-lived transitions in
# the intermediate state and can enforce dwell time in the end-state
# before a valid transition is detected. This script builds on the principles
# of the Dwell_Counter_State_Labels_VariableMem.py and
# Dwell_Counter_State_Labels.py scripts.
#
# Example: Given the input state stream where we intend to detect 1 2 3:
#          3
#          1
#          1
#          2
#          2
#          2
#          4
#          2
#          3
#          3
#          3
#          2
#          This script would still detect a 1 2 3 transition if we had
#          enforced a label_memory of 3 and a threshhold >= 0.8
#
# By Chris Ing, 2013 for Python 2.7
#
################################################################################
from Ion_Preprocessor import *
from re import match
from argparse import ArgumentParser
from collections import defaultdict
from itertools import chain, repeat, islice, groupby

# This helper function will allow me to iterate over a fixed window
# http://stackoverflow.com/q/6998245/1086154
def window(seq, size=2, fill=0, fill_left=False, fill_right=False):
    """ Returns a sliding window (of width n) over data from the iterable:
      s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...
    """
    ssize = size - 1
    it = chain(
      repeat(fill, ssize * fill_left),
      iter(seq),
      repeat(fill, ssize * fill_right))
    result = tuple(islice(it, size))
    if len(result) == size:  # `<=` if okay to return seq if len(seq) < size
        yield result
    for elem in it:
        result = result[1:] + (elem,)
        yield result

# This function will extract all transitions (including fast recrossings)
# and return a dictionary with key of type "regex_label1"-"regex_label2"
# and the value is the number of those transitions observed.
def state_transitions(data_floats, data_regex,
                      time_col=0, time_increment=1,
                      traj_col=11, verbose=False):

    transition_totals = defaultdict(int)

    # Extract the useful data out of these large lists
    data_paired = []
    for line, extra in zip(data_floats,data_regex):
        data_paired.append([line[time_col],line[traj_col]]+[extra])

    # n[0] is time, n[1] is trajectory, n[2] is regex_state
    for n, nplus1 in window(data_paired,2):
        if verbose:
            print (n, nplus1, n[1] == nplus1[1],
                  nplus1[0] - n[0] == time_increment)
        # Verify there is no missing timestep and that the run number
        # has not changed, also verify that it's a transition.
        if (n[1] == nplus1[1] and
            nplus1[0] - n[0] == time_increment and
            n[2] != nplus1[2]):
            #print N[1], Nplus1[1]
            transition_totals[str(n[2])+"-"+str(nplus1[2])] += 1

    return transition_totals

# This function processes a state stream and turns it into a list of the
# format: time traj state_id time_spent_in_that_state
# Using this format it is trivial to exclude low occupancy intermediates
# and extract 3-state permeation events. It is up to the user to interpret
# which of the 3-state events are biologically revelant and indicate a
# permeation event.
def state_intermediate_transitions(data_floats, data_regex,
                      time_col=0, time_increment=1,
                      traj_col=11, verbose=False):

    for line, extra in zip(data_floats,data_regex):
        data_paired.append([line[time_col],line[traj_col]]+[extra])

    # This time we'll group the state stream into a list of lists
    # with the magic of itertools: http://stackoverflow.com/a/7025601/1086154
    grouped_regex = [list(g) for k, g in groupby(data_regex)]

    for chunk in grouped_regex

    transition_totals = defaultdict(int)

    data_lumped = []
    data_paired = []
    for line, extra in zip(data_floats,data_regex):
        data_paired.append([line[time_col],line[traj_col]]+[extra])

    temp_state = 0
    # n[0] is time, n[1] is trajectory, n[2] is regex_state
    for n, nplus1 in window(data_paired,2):
        if n[1]


        if (n[1] != nplus1[1] and
            nplus1[0] - n[0] == time_increment and
            n[2] != nplus1[2]):

            data_lumped.append([n[1]-time_entered, n[0], n[2]])
            time_entered = nplus1[0]



        data_paired.append([line[time_col],line[traj_col]]+[extra])

if __name__ == '__main__':
    parser = ArgumentParser(
    description='This script takes regex expressions for a state label\
    and outputs how many of your states are classified by that label\
    as well as the population of those states in the dataset')

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
    '-sc', dest='sort_cut', type=float, default=0.0,
    help='a value on the sort_col range to classify zero coordinated data')
    parser.add_argument(
    '-sf', dest='sf_col', type=int, nargs="+", default=[5,6],
    help='the coordination integer columns that define the selectivity filter')
    parser.add_argument(
    '-n', dest='num_ions_map', type=int, nargs="+", default=[1,2,2,3,0],
    help='list of integer ion counts in SF for each regex value + 1 for extra')
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
    '-i', dest='regex', type=str, nargs="+", required=True,
    help='a list of regex values in quotes')

    '''
    # These are specific to this script
    parser.add_argument(
    '-b', dest='regex_begin', type=str, required=True,
    help='a regex value to start in')
    parser.add_argument(
    '-i', dest='regex_int', type=str, default="",
    help='a regex value as the intermediate state, optional')
    parser.add_argument(
    '-e', dest='regex_end', type=str, required=True,
    help='a regex value to end in')
    '''
    parser.add_argument(
    '-thres', dest='threshold', type=float, default=0.9,
    help='float from 0 to 1 that defines the occupancy % in the intermediate')
    parser.add_argument(
    '-endmem', dest='end_memory', type=int, default=2,
    help='float from 0 to 1 that defines the occupancy % in the intermediate')
    parser.add_argument(
    '-timecol', dest='time_col', type=int, default=0,
    help='a zero inclusive column number with the frame/timestep number')
    parser.add_argument(
    '-dt', dest='time_increment', type=int, default=1,
    help='the difference in the time column between steps')
    args = parser.parse_args()

    data_f, data_f_padded = process_input(filenames=args.filenames,
                                          num_cols=args.num_cols,
                                          max_ions=args.max_ions,
                                          remove_frames=args.remove_frames,
                                          traj_col=args.traj_col,
                                          sort_col=args.sort_col,
                                          add_time=args.add_time)

    data_f_regex = regex_columns(data_f_padded, regex_strings=args.regex,
                                num_cols=args.num_cols,
                                sort_col=args.sort_col,
                                sort_cut=args.sort_cut,
                                sf_col=args.sf_col,
                                max_ions=args.max_ions)

    print "State Transition Counting"
    print state_transitions(data_f_padded, data_f_regex,
                            time_col=args.time_col,
                            time_increment=args.time_increment,
                            traj_col=args.traj_col)

