#!/usr/bin/python
###############################################################################
#
# This script outputs both the rotameric state as a function of time for
# time series plotting and rotamer distributions for error bar analysis
# of rotameric states per trajectory. This script is directly analogous to
# the Occupancy Vs Time scripts. Unfortunately, these scripts are based on
# the assumption that only 2 rotamer states exist (0 and 1) and though
# the rotamer_preprocessor supports multiple cuts, this function doesn't yet
# support that in a useful way.
#
# Example: For a dihedral vs time and state-stream merged data like:
#          25000.0 2.0 ... 1.0 1.0 1.0 0.0
#          25001.0 2.0 ... 1.0 0.0 1.0 0.0
#          25002.0 2.0 ... 1.0 0.0 1.0 1.0
#
#          The following command would remove 2000 lines from the input
#          and produce a number of plots and statistical output.
#          python Rotamer_Vs_Time.py -f f1.out f2.out -t 9
#                                    -x1 0 2 3 4 -x2 1 3 4 5
#                                    -x2_cut 180
#                                    -remove 2000
#
# By Chris Ing, 2013 for Python 2.7
#
###############################################################################
from argparse import ArgumentParser
from numpy import mean
from scipy.stats import sem
from collections import defaultdict
from ChannelAnalysis.RotamerAnalysis.Preprocessor import *

# This function is useful for counting the ratio of rotameric states
# for each of the input trajectory files with respect to the last state.
# So you might get a traj_id 10 that has a state ratio of 8:1 which means
# that there was 8 times as much in state 1.
def rotamer_counter_per_res(data_lines, data_states, traj_col=11, prefix=None):

    # This is an epic datatype that I will use to quickly build a
    # dict of dicts where the 1st key is a trajectory number
    # and the second key is the ion count and the value is a
    # count of how many times that ion count was observed.
    count_totals = defaultdict(lambda: defaultdict(int))

    # This is the return list in the format:
    # traj_id d/u_ratio uncertainty
    for line, states in zip(data_lines, data_states):
        traj_id = line[traj_col]
        for state in states:
            count_totals[traj_id][state] += 1

    # D/U ratios for all trajectories
    rotamer_ratios = []
    for traj_id, count_dict in count_totals.iteritems():
        # We want the highest dunking state (1 in most cases) and to extract
        # the count for that state as the temp_max variable
        temp_max = float(sorted(count_dict.iteritems())[-1][-1])

        temp_row = [traj_id]
        for state, state_count in sorted(count_dict.iteritems()):
            temp_row.append(state_count/temp_max)
        rotamer_ratios.append(temp_row)

    # Hack to tranpose list of lists and iterate over rotamer states
    flipped_rows = zip(*rotamer_ratios)
    temp_mean = ["MEAN"]
    temp_sem = ["STDERR"]
    # The first flipped_row entry is a list of traj_id's so we skip it
    for rotamer_state in flipped_rows[1:]:
        temp_mean.append(mean(rotamer_state))
        temp_sem.append(sem(rotamer_state))
    rotamer_ratios.append(temp_mean)
    rotamer_ratios.append(temp_sem)

    return rotamer_ratios

# This function counts the number of dunking states at each step for each
# of the dihedral columns classified using the label_states function
# and separates them on a basis of trajectory id. This allows for a
# statistical measure of the distribution of dunking states much like
# the channel and selectivity filter occupancy functions.
def rotamer_counter(data_lines, data_states, traj_col=11, prefix=None):

    # This is an epic datatype that I will use to quickly build a
    # dict of dicts where the 1st key is a trajectory number
    # and the second key is the ion count and the value is a
    # count of how many times that ion count was observed.
    count_totals = defaultdict(lambda: defaultdict(int))

    for line, states in zip(data_lines, data_states):
        traj_id = line[traj_col]
        state_total = sum(states)
        count_totals[traj_id][state_total] += 1

    print count_totals
    return count_totals_to_percents(count_totals)

# This is a helper function that takes the datatype generated in
# *_counter functions (trajnum dict -> occupancy_id -> integer counts)
# and converts this to populations in a list. Num ions map is
# useful the occupancy_id's represent distinct numbers of ions
# in the selectivity filter.
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
    all_weighted_avgs=[]
    weighted_avgs_by_occid=defaultdict(list)
    for traj_id, percents in ion_count_percents.iteritems():
        temp_weighted_avg = 0
        for occ_id, percent in enumerate(percents):
            x = ion_count_indices[traj_id][occ_id]*percent
            temp_weighted_avg += x
            weighted_avgs_by_occid[occ_id].append(x)
        all_weighted_avgs.append(temp_weighted_avg)

    for occ_id, weight_avg in weighted_avgs_by_occid.iteritems():
        ion_count_percents['MEAN'].append(mean(weight_avg))
        ion_count_indices['MEAN'].append(occ_id)
        ion_count_percents['STDERR'].append(sem(weight_avg))
        ion_count_indices['STDERR'].append(occ_id)

    ion_count_percents['MEAN'].append(mean(all_weighted_avgs))
    ion_count_indices['MEAN'].append('ALL')
    ion_count_percents['STDERR'].append(sem(all_weighted_avgs))
    ion_count_indices['STDERR'].append('ALL')

    return (dict(ion_count_percents), dict(ion_count_indices))

# This writes the total of all rotamer states as a function of time
# for use of time series plotting purposes.
def compute_rotamer_vs_time(data_lines, data_states, traj_col, prefix=None):

    # This a dictionary of file streams that will be used for output
    # when prefix is assigned.
    count_files = {}

    # These are dictionaries of dictionaries where the key is a trajectory
    # number and the list is the computed occupancy count, or frame number
    rotstate_per_traj = defaultdict(list)
    time_per_traj = defaultdict(list)

    # This is an epic datatype that I will use to quickly build a
    # dict of dicts where the 1st key is a trajectory number
    # and the second key is the ion count and the value is a
    # count of how many times that ion count was observed.
    count_totals=defaultdict(lambda: defaultdict(int))

    for line, states in zip(data_lines, data_states):
        traj_id = line[traj_col]
        state_total = sum(states)
        count_totals[traj_id][state_total] += 1

        if prefix != None:
            if traj_id not in count_files:
                count_files[traj_id] = \
                             open(prefix+"_total_n"+str(traj_id),"w")
            count_files[traj_id].write(str(line[0])+" "+
                                           str(traj_id)+" "+
                                           str(state_total)+
                                           "\n")

        rotstate_per_traj[traj_id].append(state_total)
        time_per_traj[traj_id].append(float(line[0]))

    # Close filestreams.
    if prefix != None:
        for key in count_files.keys():
            count_files[key].close()

    return (dict(rotstate_per_traj), dict(time_per_traj))

# Returns the population percentages across the entire dataset from the
# label_states array.
def label_statistics(data_states):
    cols_per_line = len(data_states[0])
    counts_per_state = defaultdict(float)
    for line in data_states:
        for state in line:
            counts_per_state[state] += 1.0

    for state in counts_per_state.keys():
        counts_per_state[state] /= len(data_states)*cols_per_line

    return dict(counts_per_state)

if __name__ == '__main__':
    parser = ArgumentParser(
    description='Produces rotamer state timeseries datafiles')

    parser.add_argument(
    '-f', dest='filenames', type=str, nargs="+", required=True,
    help='a filename of coordination data from MDAnalysis trajectory data')
    parser.add_argument(
    '-x1', dest='chi1_cols', type=int, nargs="+", required=True,
    help='column numbers in the input that denote chi1 values')
    parser.add_argument(
    '-x2', dest='chi2_cols', type=int, nargs="+", required=True,
    help='column numbers in the input that denote chi2 values')
    parser.add_argument(
    '-remove', dest='remove_frames', type=int, default=0,
    help='this is a number of frames to remove from the start of the data')
    parser.add_argument(
    '-div', dest='dividers', type=float, nargs="+", default=[180],
    help='slices in angle space that label dunking states (<180 = 0, >180 = 1')
    parser.add_argument(
    '-t', dest='traj_col', type=int, default=11,
    help='a zero inclusive column number that contains the run number')
    parser.add_argument(
    '-o', dest='outfile', type=str, default=None,
    help='the file to output the sorted padding output of all input files')
    args = parser.parse_args()

    data_f_dunk = process_rotamers(filenames=args.filenames,
                                  chi1_cols=args.chi1_cols,
                                  chi2_cols=args.chi2_cols,
                                  remove_frames=args.remove_frames,
                                  traj_col=args.traj_col)

    # Same thing but now we pass the SF column list.
    print "Dunking states using the dividers", args.dividers
    data_f_states = label_states(data_f_dunk, args.chi2_cols, args.dividers)

    print "Computing dunking populations"
    print rotamer_counter(data_f_dunk, data_f_states, traj_col=args.traj_col)

    print compute_rotamer_vs_time(data_f_dunk, data_f_states, traj_col=args.traj_col,
                          prefix="rotamer")

    print "Computing dunking counts per residue"
    print rotamer_counter_per_res(data_f_dunk, data_f_states,
                                  traj_col=args.traj_col)