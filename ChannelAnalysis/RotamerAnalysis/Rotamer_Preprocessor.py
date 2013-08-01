#!/usr/bin/python
###############################################################################
#
# This script prepares multiple files of rotamer analysis output from
# MDAnalysis scripts in order to facilitate numerous Chi1-Chi2 analysis
# protocols including rotameric kinetics analysis and subsequent exponential
# fitting. This script in particular prepares a list of lists datatype
# that will be used as intput for these rotameric analysis functions.
#
# Example: For 10-column data with this type (described elsewhere):
#          24.0 170.2 63.7 -74.5 -74.4 -143.9 -88.4 -143.6 -73.7 10.0
#          25.0 157.4 78.0 -88.5 -73.3 -144.4 -73.7 -151.2 -55.1 10.0
#
#          The following command will load two datafiles, compute the
#          actual dihedral angles (shifted from MDAnalysis to 0-360),
#          label rotamer states based on a x2 cutoff and remove 2000 lines
#
#          python Rotamer_Preprocessor.py -f f1.out f2.out -t 9
#                                         -x1 0 2 3 4 -x2 1 3 4 5
#                                         -x2_cut 180
#                                         -remove 2000
#
# This script typically is the third step in a larger analysis pipeline.
# As far as script history, it's a combination of ChiToDunkCount.py
#
# By Chris Ing, 2013 for Python 2.7
#
###############################################################################
from argparse import ArgumentParser
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

# The meat of the file processing, nothing fancy though, just shifting
# negative degrees to positive in the range of 0-360. Chi2 columns
# is actually optional if your sidechain doesn't have it just omit.
def process_rotamers(filenames, chi1_cols=[], chi2_cols=[],
                    remove_frames=0, traj_col=11, time_col=0):

    # The output list
    data_floats = []

    # Parse the input file and split and float all the columns.
    # This is the required format for all the functions in this
    # file.
    for filename in filenames:
        with open(filename,"r") as data_file:
            data_raw = data_file.readlines()[remove_frames:]
            data_raw_split = [line.strip().split() for line in data_raw]

            for line in data_raw_split:
                # Exclude a header proceeded by an exclamation, another
                # weird convention secretly adopted...
                if line[0] != "!":
                    time_val = int(float(line[time_col]))
                    traj_val = int(float(line[traj_col]))

                    chi_shifted = []
                    chi_sorted_cols = sorted(chi1_cols+chi2_cols)
                    for col in chi_sorted_cols:
                        raw_chi = -1.0*float(line[col])
                        if raw_chi <= 0:
                            chi_shifted.append(raw_chi + 360.0)
                        else:
                            chi_shifted.append(raw_chi)

                    # I had this weird convention where the traj_num
                    # was at the end of the rows in the input file, but now
                    # it's at the start. Enragingly so. Legacy support.
                    if (all([traj_col > x_col for x_col in chi_sorted_cols])):
                        data_floats.append([time_val]+
                                           chi_shifted +
                                           [traj_val])
                    else:
                        data_floats.append([time_val, traj_val]+
                                           chi_shifted)

    return data_floats

# Given processed dunking data, this uses state dividers in whatever
# data range is given in the dihedral_columns and classifies it into integer
# states. Example: given the divider 60 and 300, anything lower than 60
# would be state 0, between 60 and 300 would be state 1, and state 2 for > 300.
def label_states(data_floats, dihedral_cols, dividers=[180]):

    # This is the highest or lowest value in the dihedral_cols range (+/-)
    edge = 9999

    data_output = []
    for line in data_floats:
        dihedral_vals = [line[col] for col in dihedral_cols]
        dihedral_states = []
        for dihedral in dihedral_vals:
            # Using a rolling window helper function to detect where you are
            # with an arbitrary number of state_dividers is a bit overkill
            # but it works.
            for id, range_vals in enumerate(window([-edge]+dividers+[edge])):
                if range_vals[0] < dihedral <= range_vals[1]:
                    dihedral_states.append(id)

        # If this assert fails, you have data outside the minmax range.
        assert(len(dihedral_states) == len(dihedral_vals))
        data_output.append(dihedral_states)

    return data_output

if __name__ == '__main__':
    parser = ArgumentParser(
    description='This script parses dihedral angles for multiple subunits, \
    preparing input for subsequent rotameric analysis.')

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
