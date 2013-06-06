#!/usr/bin/python

################################################################################
#
# This script parses a state stream and detects transitions that
# satisfy a particular regular expression for the state in start,
# N intermediates, and an end state. The script can exclude short-lived
# states and function with multiple trajectory datafiles as input.
# This script builds on the principles of the
# Dwell_Counter_State_Labels_VariableMem.py and Dwell_Counter_State_Labels.py
# scripts.
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
# and return a list with key of type "regex_label1"-"regex_label2"
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

    return transition_totals.items()

# This function processes a state stream and interally turns it into a
# list of the format: time traj state_id time_spent_in_that_state
#
# Using this format it is trivial to exclude low occupancy intermediates
# and extract N-state permeation events. It is up to the user to interpret
# which of the N-state events are biologically revelant and indicate a
# permeation event. The intermediates parameter adjusts the N variable
# although typically 1 is sufficient.
# The datatype returned is a list of tuples where the first tuple index
# is the transition id given by a sequence of hypen separated regex ids
# and the value is the transition count.
def state_intermediate_transitions(data_floats, data_regex,
                      time_col=0, time_increment=1, intermediates=1,
                      traj_col=11, cutoff=1, verbose=False):

    # Extract the useful data out of these large lists with the traj
    # as the key and the list as the timesteps.
    data_paired_by_traj = defaultdict(list)
    for line, extra in zip(data_floats,data_regex):
        traj = line[traj_col]
        time = line[time_col]
        data_paired_by_traj[traj].append([time,traj,extra])

    data_grouped = []
    # This time we'll group the state stream into a list of lists grouped by
    # regex_id with the magic of itertools.groupby:
    # http://docs.python.org/2/library/itertools.html#itertools.groupby
    for traj, timesteps in data_paired_by_traj.iteritems():
        for key, group in groupby(timesteps, key=lambda col: col[2]):
            group_list = list(group)
            group_len = len(group_list)*time_increment

            # Here we can impose a cutoff that removes low-population
            # states in between large intermediates.
            if group_len > cutoff:
                data_grouped.append(group_list[0]+[group_len])
            #print group_list[0]+[group_len]

    # Any cutoff value other than 0 will create gaps in the stream
    # and these need to be corrected. The following block of code
    # will look for gaps and adjust state lifetimes accordingly.
    data_collapsed = []
    # Make sure we have at least 1 group and then extract it's regex ID.
    assert len(data_grouped) > 1
    previous_group = data_grouped[0]
    start_group = data_grouped[0]
    temp_dwell = 0
    for current_group in data_grouped:
        # If there is a change of traj num since the past step
        # then add the final dwell_time value before resetting
        if previous_group[1] != current_group[1]:
            temp_dwell += previous_group[3]
            # We only take the first 3 columns because we're rewriting
            # the last column.
            data_collapsed.append(start_group[:3]+[temp_dwell])
            temp_dwell = 0
            start_group = current_group
        # If there is a change of id since the past step
        # use the difference in time values to correct the dwell time
        elif previous_group[2] != current_group[2]:
            temp_dwell += current_group[0]-previous_group[0]
            data_collapsed.append(start_group[:3]+[temp_dwell])
            temp_dwell = 0
            start_group = current_group
        else:
            temp_dwell += current_group[0]-previous_group[0]
            #print current_group[0], temp_dwell
        previous_group = current_group

    temp_dwell += previous_group[0]-start_group[0]
    data_collapsed.append(start_group[:3]+[temp_dwell])

    # This dictionary is returned by the function and summarizes the observed
    # transitions.
    transition_totals = defaultdict(int)
    # intermediates+2 is used because we have start and end states to consider
    for group_window in window(data_collapsed,intermediates+2):
        # n[0] is time, n[1] is traj_id, n[2] is regex_id, n[3] is dwell_time
        group_regex_ids = [str(temp_group[2]) for temp_group in group_window]
        # Confirm that a real transition ocurred oppossed to a back-crossing
        if len(group_regex_ids) == len(set(group_regex_ids)):
            transition_totals["-".join(group_regex_ids)] += 1

    return transition_totals.items()

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

    print "State Transition Counting"
    print state_transitions(data_f_padded, data_f_regex,
                            time_col=args.time_col,
                            time_increment=args.time_increment,
                            traj_col=args.traj_col)

    print "State Transition w/ Intermediate Counting"
    print state_intermediate_transitions(data_f_padded, data_f_regex,
                                         time_col=args.time_col,
                                         time_increment=args.time_increment,
                                         traj_col=args.traj_col)
