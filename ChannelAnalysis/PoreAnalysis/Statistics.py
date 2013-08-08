#!/usr/bin/python
###############################################################################
#
# This script extracts RMSD fluctuations from chain sliding data and computes
# it on a per trajectory basis as well as over the complete dataset.
#
# Example: For 2+12-column data like (with reduced decimals for example):
#          1.0 19.0 -0.1 -0.1 0.5 -0.1 -0.1 0.5 -0.1 -0.1 0.5 0.2 0.4 -0.5
#
#          This script would output the mean deviation with standard error of
#          mean in a format like this:
#
# By Chris Ing, 2013 for Python 2.7
#
###############################################################################
from argparse import ArgumentParser
from ChannelAnalysis.PoreAnalysis import *
from numpy import mean, sqrt, array, square
from scipy.stats import sem

#a great helper function to iterate over chunks of a list
def chunker(seq, size):
    return (seq[pos:pos + size] for pos in xrange(0, len(seq), size))

# This computes the root mean square deviation of each of the sort_col
# columns of each line.
def rmsd_counter(data_lines, colskip=2, num_cols=3,
                 sort_col=2, traj_col=1, chain_num=None, prefix=None):

    rmsd_totals = defaultdict(list)
    rmsd_alltotals = []
    rmsd_means = defaultdict(int)
    rmsd_stderrs = defaultdict(int)

    # First determine the mean displacement for the entire dataset.
    traj_means = 0.0
    for line in data_lines:
        col_blocks = list(chunker(line[colskip:],num_cols))
        traj_means += mean([block[sort_col] for block in col_blocks])
    traj_means /= len(data_lines)

    for line in data_lines:
        traj_id = line[traj_col]
        temp_deviation = []

        # Split the line into chunks of size equal to num_cols
        col_blocks = list(chunker(line[colskip:],num_cols))

        # This script can calculate statistics for 1 chain, or all chains
        # averaged, depending on if the chain_num argument is set.
        if chain_num is None:
            for ion in col_blocks:
                temp_deviation.append(float(ion[sort_col]))
        else:
            temp_deviation.append(float(col_blocks[chain_num][sort_col]))
            #print float(col_blocks[chain_num][sort_col]),

        shifted_deviation = square(array(temp_deviation)-traj_means)
        # If chain_num is set, return deviation without dviding.
        if chain_num is None:
            rmsd = sqrt(sum(shifted_deviation))/len(col_blocks)
        else:
            rmsd = sqrt(shifted_deviation)

        rmsd_totals[traj_id].append(rmsd)

    for traj_id in rmsd_totals.keys():
        rmsd_means[traj_id] = mean(rmsd_totals[traj_id])
        rmsd_stderrs[traj_id] = sem(rmsd_totals[traj_id])
        rmsd_alltotals.append(mean(rmsd_totals[traj_id]))

    rmsd_means["ALL"] = mean(rmsd_alltotals)
    rmsd_stderrs["ALL"] = sem(rmsd_alltotals)

    return (dict(rmsd_means), dict(rmsd_stderrs))

if __name__ == '__main__':
    parser = ArgumentParser(
    description='This script computes statistics on pore sliding')

    parser.add_argument(
    '-f', dest='filenames', type=str, nargs="+", required=True,
    help='a filename of pore sliding data from MDAnalysis')
    parser.add_argument(
    '-c', dest='num_cols', type=int, default=3,
    help='the number of columns per ion in the input, typically x,y,z=3')
    parser.add_argument(
    '-remove', dest='remove_frames', type=int, default=0,
    help='this is a number of frames to remove from the start of the data')
    parser.add_argument(
    '-s', dest='sort_col', type=int, default=2,
    help='a zero inclusive column number to sort your row on, typically z=2')
    parser.add_argument(
    '-t', dest='traj_col', type=int, default=1,
    help='a zero inclusive column number that contains the run number')
    args = parser.parse_args()

    sf_processed = process_channelatoms(args.filenames,
                                        remove_frames=args.remove_frames)

    if True:
        print "All Chains Mean Deviation"
        sliding_stats = rmsd_counter(sf_processed,
                           num_cols=args.num_cols,
                           traj_col=args.traj_col,
                           sort_col=args.sort_col)

        for traj_id in sliding_stats[0].keys():
            print traj_id,
            print sliding_stats[0][traj_id], " +- ",
            print sliding_stats[1][traj_id]

    else:
        print "Chain A Mean Deviations"
        sliding_stats = rmsd_counter(sf_processed,
                           num_cols=args.num_cols,
                           traj_col=args.traj_col,
                           sort_col=args.sort_col,
                           chain_num=0)

        for traj_id in sliding_stats[0].keys():
            print traj_id,
            print sliding_stats[0][traj_id], " +- ",
            print sliding_stats[1][traj_id]

        print "Chain B Mean Deviations"
        sliding_stats = rmsd_counter(sf_processed,
                           num_cols=args.num_cols,
                           traj_col=args.traj_col,
                           sort_col=args.sort_col,
                           chain_num=1)

        for traj_id in sliding_stats[0].keys():
            print traj_id,
            print sliding_stats[0][traj_id], " +- ",
            print sliding_stats[1][traj_id]

        print "Chain C Mean Deviations"
        sliding_stats = rmsd_counter(sf_processed,
                           num_cols=args.num_cols,
                           traj_col=args.traj_col,
                           sort_col=args.sort_col,
                           chain_num=2)

        for traj_id in sliding_stats[0].keys():
            print traj_id,
            print sliding_stats[0][traj_id], " +- ",
            print sliding_stats[1][traj_id]

        print "Chain D Mean Deviations"
        sliding_stats = rmsd_counter(sf_processed,
                           num_cols=args.num_cols,
                           traj_col=args.traj_col,
                           sort_col=args.sort_col,
                           chain_num=3)

        for traj_id in sliding_stats[0].keys():
            print traj_id,
            print sliding_stats[0][traj_id], " +- ",
            print sliding_stats[1][traj_id]
