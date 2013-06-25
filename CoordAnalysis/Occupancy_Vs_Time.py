#!/usr/bin/python
###############################################################################
#
# This script performs time-series data on one, multiple, or merged
# data from the Row_Preprocessor.py script. The purpose of this tool
# is to extract several important properties of this data as a function
# of time, as well as some basic statistical properties over
# the given ensembles. In terms of script history, this combines
# ion_occupancy_vs_time.py, Compute_Average_Sod_in_SF.py
#
# Example: For 13/26/39/...-column data with type like:
#          1.0 -0.13 -0.193 0.522 0.0 0.0 0.0 0.0 0.0 2.0 9.0 2.0 1748.0
#
#          The following command would remove 200 lines from the input
#          and compute the SF occupancy using columns 5,6 of each ion
#          grouping:
#          python Occupancy_Vs_Time.py -f f1.out -m 3 -c 13 -remove 200 -sf 5 6
#
# By Chris Ing, 2013 for Python 2.7
#
###############################################################################
from argparse import ArgumentParser
from numpy import mean
from scipy.stats import sem
from collections import defaultdict
from Ion_Preprocessor import *

# This counts of the number of ions at each step and bins these based on
# integer occupancy values. This function utilizes the traj_col variable
# assuming that multiple raw data files are passed as input and returns
# the mean ion occupancy value across all trajectories as well as the
# standard error of mean. If a prefix is specified, the entire dataset
# is split into datafiles of different ion occupancy.
# If the sf_col list is set, it will get for coordination of ions
# in any of those columns and only tabulate occupancy based on satisfying
# that criteria first.
def occ_counter(data_lines, num_cols=13, traj_col=11,
                pad_col=4, sf_col=[5,6], prefix=None):

    # This is an epic datatype that I will use to quickly build a
    # dict of dicts where the 1st key is a trajectory number
    # and the second key is the ion count and the value is a
    # count of how many times that ion count was observed.
    count_totals=defaultdict(lambda: defaultdict(int))

    # This a dictionary of file streams that will be used for output
    # when prefix is assigned.
    count_files={}

    for line in data_lines:
        traj_id = line[traj_col]
        temp_ion_count = 0
        for ion in chunker(line,num_cols):
            # In the case that the ion grouping has a hypen
            # we know that's a padded column and must be excluded.
            if ion[pad_col] != "-":
                if sf_col != []:
                    if any([ion[sf_colid] > 0.0 for sf_colid in sf_col]):
                        temp_ion_count += 1
                else:
                    temp_ion_count += 1

        # Good old defaultdict, no need to initialize!
        count_totals[traj_id][temp_ion_count] += 1

        # File streams for splitting purposes.
        if prefix != None:
            if temp_ion_count not in count_files:
                count_files[temp_ion_count] = \
                             open(prefix+"_split"+str(temp_ion_count),"w")
            else:
                count_files[temp_ion_count].write(
                    " ".join([str(col) for col in line]))
                count_files[temp_ion_count].write("\n")

    if prefix != None:
        for key in count_files.keys():
            count_files[key].close()

    # Return the list of list, the mean and standard error of mean
    # for each trajectory in the input.
    return count_totals_to_percents(count_totals)

# This is a helper function that takes the datatype generated in
# *_counter functions (trajnum dict -> occupancy_id -> integer counts)
# and converts this to populations in a list. Num ions map is
# useful the occupancy_id's represent distinct numbers of ions
# in the selectivity filter.
def count_totals_to_percents(count_totals, num_ions_map=[]):

    # Very ugly line that simply finds the maximum key for ion counts above
    max_ions = max([traj for count_dict in count_totals.values()
                    for traj in count_dict])

    if not num_ions_map:
        num_ions_map = range(max_ions+1)
    assert len(num_ions_map) == max_ions+1, "Insufficient number of elements"

    # This list has elements equal to the number of timesteps
    # where each sublist contains max_ions rows and percentages
    # under each column for that occupancy.
    traj_total_percent = []

    # Yet another average, here we want to know the mean and stdev occupancy
    # of 1-ion, 2-ion, across all trajectories. They key is occupancy number
    # and the value is a list of occupancies in percent.
    traj_occ_averages = defaultdict(list)

    for traj_id, count_dict in count_totals.iteritems():
        traj_total_lines = float(sum(count_dict.values()))

        temp_line = []
        temp_line.append(traj_id)
        for traj_index in range(max_ions+1):
            temp_percent = count_dict[traj_index]/traj_total_lines
            temp_line.append(temp_percent)
            traj_occ_averages[traj_index].append(temp_percent)

        # Finally we append the last column which is the ion occupancy
        # average for each trajectory
        temp_average=0.
        for occupancy,ion_count in zip(temp_line[1:],num_ions_map):
            temp_average += occupancy*ion_count

        traj_total_percent.append(temp_line+[temp_average])

        # We also store the weighted average in the same dictionary
        # Note that max_ions+1 index doesn't exist in the loop
        # above, we're using it to store this new value.
        traj_occ_averages[max_ions+1].append(temp_average)

    mean_occupancy =  [mean(traj_occ_averages[traj_index]) \
                       for traj_index in range(max_ions+2)]

    stderr_occupancy = [sem(traj_occ_averages[traj_index]) \
                        for traj_index in range(max_ions+2)]

    traj_total_percent.append(["MEAN"]+mean_occupancy)
    traj_total_percent.append(["STDERR"]+stderr_occupancy)

    return traj_total_percent

# Similar to the function above, but this script outputs
# number of ions as a function of time to output files.
def compute_occ_vs_time(data_lines, num_cols=13,
                      traj_col=11, pad_col=4, sf_col=[], prefix=None):

    # This a dictionary of file streams that will be used for output
    # when prefix is assigned.
    count_files={}

    # This is a dictionary of dictionaries where the key is a trajectory
    # number and the list is the computed occupancy count
    occ_per_traj_vs_time=defaultdict(list)

    for line in data_lines:
        traj_id = line[traj_col]
        temp_ion_count = 0
        for ion in chunker(line,num_cols):
            # In the case that the ion grouping has a hypen
            # we know that's a padded column and must be excluded.
            if ion[pad_col] != "-":
                if sf_col != []:
                    if any([ion[sf_colid] > 0.0 for sf_colid in sf_col]):
                        temp_ion_count += 1
                else:
                    temp_ion_count += 1

            if prefix != None:
                if traj_id not in count_files:
                    count_files[traj_id] = \
                                 open(prefix+"_occupancy_n"+str(traj_id),"w")
                count_files[traj_id].write(str(ion[0])+" "+
                                           str(traj_id)+" "+
                                           str(temp_ion_count)+
                                           "\n")

            occ_per_traj_vs_time[traj_id].append(temp_ion_count)

    # Close filestreams.
    if prefix != None:
        for key in count_files.keys():
            count_files[key].close()

    return occ_per_traj_vs_time.items()


if __name__ == '__main__':
    parser = ArgumentParser(
    description='This script computes channel and sf occupancy values\
    and statistics for time series coordination input')

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
    '-sf', dest='sf_col', type=int, nargs="+", default=[5,6],
    help='the coordination integer columns that define the selectivity filter')
    parser.add_argument(
    '-o', dest='outfile', type=str, default=None,
    help='the file to output the sorted padding output of all input files')
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

    # Running ion_counter on the original data_float list will give
    # an honest estimate of the ionic population in the channel.
    # This can also be run after you've done the padding step (which
    # will omit ions > max_ions, fixing your number of columns in the
    # datafile)
    print "Channel Occupancy"
    print occ_counter(data_f, num_cols=args.num_cols,
                      traj_col=args.traj_col, sf_col=[], prefix="chanocc")

    print compute_occ_vs_time(data_f, num_cols=args.num_cols,
                      prefix="chanocc")

    # Same thing but now we pass the SF column list.
    print "SF Occupancy using columns", args.sf_col
    print occ_counter(data_f, num_cols=args.num_cols,
                      traj_col=args.traj_col, sf_col=args.sf_col,
                      prefix="sfocc")

    print compute_occ_vs_time(data_f, num_cols=args.num_cols,
                          sf_col=args.sf_col, prefix="sfocc")
