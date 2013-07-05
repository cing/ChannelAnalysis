#!/usr/bin/python
###############################################################################
#
# This script produces state labels for all ionic states in MDAnalysis output
# based on coordination integers and groups these according to a passed
# regular expression. This script is really just a glorified rewrite of
# Regular_Expression_Checker_w_Z.py and Regular_Expression_and_Label_Printer.py
#
# Example: For 13/26/39/...-column data with type like:
#          1.0 -0.13 -0.193 0.522 0.0 1.0 0.0 0.0 0.0 2.0 9.0 2.0 1748.0
#
#          python Regex_Grouping.py -f f1.out f2.out -m 3 -c 13 -remove 2000
#                                  -i "(.[^-0+]|[^-0+].)[-0][-0][-0][-0]"
#                                  "\+\+[-0][-0][-0][-0]"
#                                  "(.[^-+0]|[^-+0].)(.[^-+0]|[^-+0].)[-0][-0]"
#                                  "\+\+(.[^-0+]|[^-0+].)[-0][-0]"
#                                  -n 1 0 2 1 0 -sf 5 6
#
#          This would read in the finals columns 5 and 6, produce state labels
#          like ++0100, 101010, ++00--, 1010-- and then match them to the
#          four passed regular expressions. These regular expressions have
#          SF occupancy specified by the -n argument. This would produce
#          data like this where there is N+2 lists where N is the number of
#          files passed in the -f argument above:
#
#          [[2.0, 0.43, 0.53, 0.01, 0.02, 0.00, 0.48],
#           [3.0, 0.13, 0.29, 0.16, 0.40, 0.00, 0.87],
#           ['MEAN', 0.28, 0.41, 0.09, 0.21, 0.00, 0.67],
#           ['STDERR', 0.14, 0.11, 0.07, 0.18, 0.00, 0.19]]
#
# By Chris Ing, 2013 for Python 2.7
#
###############################################################################
from argparse import ArgumentParser
from numpy import mean
from scipy.stats import sem
from collections import defaultdict
from re import match
from Ion_Preprocessor import *

# This function counts regex occupancy in each trajectory of the input
# data for each of the passes regex strings. It will return the
# regex populations for each trajectory and then the associated stats.
def regex_counter(data_floats, data_regex, traj_col=11):

    # This is an epic datatype that I will use to quickly build a
    # dict of dicts where the 1st key is a trajectory number
    # and the second key is the regex index and the value is a
    # count of how many times that ion count was observed.
    count_totals=defaultdict(lambda: defaultdict(int))

    # the data_regex datatype is a list of tuples: (state_label, regex_int)
    data_regex_ids = [data[1] for data in data_regex]

    for line, regex_state in zip(data_floats,data_regex_ids):
        traj_id = line[traj_col]
        count_totals[traj_id][regex_state] += 1

    # Return the list of list, the mean and standard error of mean
    # for each trajectory in the input.
    return count_totals_to_percents(count_totals)

# This is a helper function that takes the datatype generated in
# *_counter functions (trajnum dict -> regex_id -> integer counts)
# and converts this to populations in a list without weighting like
# the occupancy count function.
def count_totals_to_percents(count_totals):

    # Here's the return datatype that stores the percentage of occupancy
    # in a given channel/sf state which can be paired with the indices
    ion_count_percents = defaultdict(list)
    ion_count_indices = defaultdict(list)
    for traj_id, count_dict in count_totals.iteritems():
        traj_total_lines = float(sum(count_dict.values()))
        for ion_state, ion_count in count_dict.iteritems():
            ion_count_percents[traj_id].append(ion_count/traj_total_lines)
            ion_count_indices[traj_id].append(ion_state)

    # Append a little statistics, sorry if this is confusing...
    avgs_by_regex=defaultdict(list)
    for traj_id, percents in ion_count_percents.iteritems():
        regex_ids = ion_count_indices[traj_id]
        for regex_id, percent in zip(regex_ids, percents):
            avgs_by_regex[regex_id].append(percent)

    for regex_id, avg in avgs_by_regex.iteritems():
        ion_count_percents['MEAN'].append(mean(avg))
        ion_count_indices['MEAN'].append(regex_id)
        ion_count_percents['STDERR'].append(sem(avg))
        ion_count_indices['STDERR'].append(regex_id)

    return (dict(ion_count_percents), dict(ion_count_indices))

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

    print "Regex Macrostate Occupancy"
    print regex_counter(data_f_padded, data_f_regex,
                        traj_col=args.traj_col)
