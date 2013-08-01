#!/usr/bin/python
###############################################################################
#
# Prepares histograms for individual rings of channel atoms based on a user-
# -defined column of the channel atom datafiles.
#
# Example: For 14-column data with this type (described elsewhere):
#
# 1 7.0 0.413 0.373 0.294 0.300 0.282 0.425 0.358 0.246 0.422 0.305 0.392 0.350
# 2 7.0 0.412 0.337 0.280 0.388 0.292 0.419 0.384 0.233 0.469 0.287 0.389 0.301
#
#          The following command will load that datafile into memory, strip
#          the first 2000 lines and produce a series of histogram datafiles
#          or return data that could be plotted accordingly using matplotlib.
#
#          python ChannelAtom_Histograms.py -f nav.n7.thr -r 4 -remove 2000
#
# By Chris Ing, 2013 for Python 2.7
#
###############################################################################
from argparse import ArgumentParser
from collections import defaultdict
from numpy import histogram, convolve, ones
from ChannelAtom_Preprocessor import *

# a great helper function to iterate over chunks of a list
def chunker(seq, size):
    return (seq[pos:pos + size] for pos in xrange(0, len(seq), size))

# a helper method for extracting a timeseries window.
def window(size):
    return ones(size)/float(size)

# This returns the sort_column as a time series, useful
# for making scatterplot time series of channel atom positions.
def compute_atom_timeseries(data_lines, sort_col, traj_col,
                            col_skip=2, num_cols=3, window_size=100):

    # These are dictionaries of dict where the key is the traj_number
    # and the subdict is ion_number and te value is a LIST of ion positions,
    # or associated time values in the case of the associated time_per_traj
    atom_pos_per_traj = defaultdict(dict)
    time_per_traj = defaultdict(dict)

    for line in data_lines:
        traj_id = line[traj_col]
        for atom_num, atom in enumerate(list(chunker(line[col_skip:],
                                                     num_cols))):
            sort_val = atom[sort_col]
            if atom_num not in atom_pos_per_traj[traj_id]:
                atom_pos_per_traj[traj_id][atom_num] = [sort_val]
                time_per_traj[traj_id][atom_num] = [line[0]]
            else:
                atom_pos_per_traj[traj_id][atom_num].append(sort_val)
                time_per_traj[traj_id][atom_num].append(line[0])

    if convolve != None:
        for t_id, atoms in atom_pos_per_traj.iteritems():
            for a_id, atom_ts in atoms.iteritems():
                atom_pos_per_traj[t_id][a_id] = list(convolve(atom_ts,
                                                     window(window_size),
                                                     'same'))

    return (dict(atom_pos_per_traj), dict(time_per_traj))

# Not a complicated function for getting histogrammed data for the sort_col
# in a particular group of data_lines. This does not distinguish between
# any of the residues in the ring, i.e. if one is protonated this will
# be lumped in all together.
def compute_allatom_histogram(data_lines, sort_col,
                             num_cols=3,
                             histmin=-1.50, histmax=1.5,
                             histbins=300, col_skip=2,
                             normed=True, prefix=None):

    # Power datatypes son. The structure is: traj_id -> ion_num -> z_vals
    atom_sortvals = []

    # These are dictionaries of lists where the key is a coord_col
    # and the list is a axial probability or associated z value.
    coord_hist_per_atom = defaultdict(list)
    z_per_atom = defaultdict(list)

    for line in data_lines:
        for atom in chunker(line[col_skip:],num_cols):
            sort_val = atom[sort_col]
            atom_sortvals.append(sort_val)

    histo, edges = histogram(atom_sortvals, range=[histmin, histmax],
                            bins=histbins, normed=normed)

    if prefix != None:
        with open(prefix+"_allatom","w") as out:
            for xval, yval in zip(edges,histo):
                out.write(str(xval)+" "+str(yval)+"\n")

    coord_hist_per_atom["ALL"].extend(list(histo))
    z_per_atom["ALL"].extend(list(edges))

    return (dict(coord_hist_per_atom), dict(z_per_atom))

if __name__ == '__main__':
    parser = ArgumentParser(
    description='This script parses input columnular ASCII data\
    of channel atoms and makes it nice and pretty for subsequent analysis.')

    parser.add_argument(
    '-f', dest='filenames', type=str, nargs="+", required=True,
    help='a filename of atom data from MDAnalysis trajectory data')
    parser.add_argument(
    '-c', dest='num_cols', type=int, default=3,
    help='the number of columns per channel atom in the input')
    parser.add_argument(
    '-cs', dest='col_skip', type=int, default=2,
    help='the number of columns per line in input that are not atom data')
    parser.add_argument(
    '-s', dest='sort_col', type=int, default=2,
    help='a zero inclusive column number to pull from each res, typically z')
    parser.add_argument(
    '-remove', dest='remove_frames', type=int, default=0,
    help='this is a number of frames to remove from the start of the data')
    args = parser.parse_args()

    data_f_processed = process_channelatoms(filenames=args.filenames,
                                            remove_frames=args.remove_frames)

    allatom_histo = compute_allatom_histogram(data_f_processed,
                                              args.sort_col,
                                              col_skip=args.col_skip,
                                              num_cols=args.num_cols)

    print allatom_histo