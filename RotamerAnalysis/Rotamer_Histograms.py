#!/usr/bin/python
###############################################################################
#
# This script generates dihedral populations of passed columns in the
# rotamer datafile. It also parses the survival probabilities of dunking states
# and outputs histograms of dunking survival. The name is kind of a lie
# since it makes 2D histograms as well, but I didn't feel like changing my
# naming convention.
#
# Example: For 13/26/39/...-column data with type like:
#          1.0 -0.13 -0.193 0.522 0.0 1.0 0.0 0.0 0.0 2.0 9.0 2.0 1748.0
#
#          The following command would remove 2000 lines from the input
#          and produce a large number of 1D Histograms to be plotted
#          in an external plotting program:
#          python Coordination_Histograms.py -f f1.out -m 3 -c 13
#                                            -remove 2000
#                                            -coord 4 5 6 7 8 9 10
#
# By Chris Ing, 2013 for Python 2.7
#
###############################################################################
from argparse import ArgumentParser
from numpy import histogram, histogram2d, log
from Rotamer_Preprocessor import *

# This writes chi1 or chi2 population densities in 1D
# (kind of a weaker version of write_rotamer_rotamer_histogram below)
# Note that it's unnormalized, but that's easy to change.
def write_rotamer_histogram(data_lines, dihedral_cols,
                             histmin=0, histmax=360, histbins=250,
                             prefix="1dhisto_chi2"):
    dihedral_vals = []
    for line in data_lines:
        dihedral_vals.extend([line[col] for col in dihedral_cols])

    histo, edges = histogram(dihedral_vals, range=[histmin, histmax],
                                     bins=histbins, normed=False)

    with open(prefix+"_dihedrals","w") as out:
        for xval, yval in zip(edges,histo):
            out.write(str(xval)+" "+str(yval)+"\n")

    return True

# This writes Chi1 vs Chi2 distributions (or vice versa) in the form of a
# 2D histogram.
def write_rotamer_rotamer_histogram(data_lines,
                                    dihedral_cols_x, dihedral_cols_y,
                                    histmin=0, histmax=360, histbins=250,
                                    prefix="2dhisto", kBT=0.596):

    dihedral_vals_x = []
    dihedral_vals_y = []
    for line in data_lines:
        dihedral_vals_x.extend([line[col] for col in dihedral_cols_x])
        dihedral_vals_y.extend([line[col] for col in dihedral_cols_y])

    histo, xedges, yedges = histogram2d(dihedral_vals_x,
                                        dihedral_vals_y,
                                        range=[[histmin,histmax],
                                               [histmin,histmax]],
                                        bins=[histbins,histbins], normed=True)

    with open(prefix+"_dihedrals_pmf","w") as out:
        for xbin in range(histbins):
            for ybin in range(histbins):
                if histo[xbin][ybin] > 0:
                    out.write(str(xedges[xbin])+" "+str(yedges[ybin])+" "+
                              str(histo[xbin][ybin])+" "+
                              str(-kBT*log(histo[xbin][ybin]))+"\n")
                else:
                    out.write(str(xedges[xbin])+" "+str(yedges[ybin])+" "+
                              str(histo[xbin][ybin])+" "+"11.0\n")
            out.write("\n")
    return True



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

    data_f_dunk = process_dunking(filenames=args.filenames,
                                  chi1_cols=args.chi1_cols,
                                  chi2_cols=args.chi2_cols,
                                  remove_frames=args.remove_frames,
                                  traj_col=args.traj_col)

    print "Writing 1D histograms chi2 population"
    write_rotamer_histogram(data_f_dunk, args.chi2_cols)

    print "Writing 2D histograms chi2 vs chi1 population"
    write_rotamer_rotamer_histogram(data_f_dunk, args.chi1_cols, args.chi2_cols)
