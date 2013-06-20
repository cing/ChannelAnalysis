#!/usr/bin/python
###############################################################################
#
# This script generates survival probabilities of rotameric states generated
# from the label_states output of Rotamer_Preprocessor. This gives you
# insight into the kinetics of rotameric motion. Survival probabilities
# are fitted to exponential functions that can be used to extract motion
# timescales.
#
# Example: For a dihedral vs time and state-stream merged data like:
#          25000.0 2.0 ... 1.0 1.0 1.0 0.0
#          25001.0 2.0 ... 1.0 0.0 1.0 0.0
#          25002.0 2.0 ... 1.0 0.0 1.0 1.0
#
#          The following command would remove 2000 lines from the input
#          and produce a number of plots and statistical output.
#          python Rotamer_Survival.py -f f1.out f2.out -t 9
#                                     -x1 0 2 3 4 -x2 1 3 4 5
#                                     -x2_cut 180
#                                     -remove 2000
#
# By Chris Ing, 2013 for Python 2.7
#
###############################################################################
from argparse import ArgumentParser
from scipy.optimize import curve_fit
from collections import defaultdict
from itertools import groupby
from Rotamer_Preprocessor import *
import numpy as np

# These are our fitting functions for multiple levels of exponential decay.
def triple_exp_func(x, a, b, c, d, e, f):
    return a * np.exp(-b * x) + c * np.exp(-d * x) + e * np.exp(-f * x)
def double_exp_func(x, a, b, c, d):
    return a * np.exp(-b * x) + c * np.exp(-d * x)
def single_exp_func(x, a, b):
    return a * np.exp(-b * x)

# This is a helper function that prepares an internal datatype for processing
# the rotamer state stream. It is used by the write_rotamer_state_survival
# function.
def process_state_lifetimes(data_lines, data_states, traj_col=11):

    # This datatype will store the state lifetime distribution with the key
    # as the rotamer id and the list of integer state lifetimes.
    data_grouped = defaultdict(list)

    # First separate the input data by the traj column, we don't want
    # to accidentally say a state is populated across two separate
    # trajectories just because they were concatted.
    data_paired_by_traj = defaultdict(list)
    for line, states in zip(data_lines,data_states):
        traj = line[traj_col]
        data_paired_by_traj[traj].append(states)

    for traj, traj_states in data_paired_by_traj.iteritems():
        # Hack to tranpose list of lists and iterate over each residue's states
        flipped_rows = zip(*traj_states)
        for rotamer, rotamer_states in enumerate(flipped_rows):
            for state_id, group in groupby(rotamer_states):
                group_list = list(group)
                group_len = len(group_list)
                data_grouped[state_id].append(group_len)

    return data_grouped

# This function creates grouped blocks of states much like Regex_Transitions
# but instead of just printing counts/rates, it builds a histogram of state
# survival times where at t=0 there is a probability of state survival of 100%
# and there is a gradual decay from there at subsequent timesteps.
# time_output_conv is a conversion from the units of the input to the output
# units (often nanoseconds)
def write_rotamer_state_survival(data_lines, data_states, survival_cut=600,
                                 traj_col=11, time_output_conv=0.02,
                                 prefix="1dhisto_rotamer_survival"):

    data_grouped = process_state_lifetimes(data_lines, data_states, traj_col)

    fit_values = []
    for state_id, lifetimes in data_grouped.items():
        # Here we make a quick 1D histogram with a bin for every
        # discrete value. Since it's a survival probability, we
        # increment all values up to t=lifetime and normalize after
        # such that the value at t=0 is 1.
        histo = [0]*(survival_cut)
        for life in lifetimes:
            # Since we're normalizing, any values greater than survival
            # cut can be excluded since we're just adding +1 to all bins.
            if life < survival_cut:
                for sublife in range(life):
                    histo[sublife] += 1

        histo = [population/float(histo[0]) for population in histo]

        # For each of the states, print the output and normalize
        # that badboy.
        with open(prefix+"_state"+str(state_id),"w") as out:
            for survival_time, population in enumerate(histo):
                out.write(str(survival_time*time_output_conv)+" "
                          +str(population)+"\n")

        # These are input values that assume a 0.1, 1, 10 nanosecond timescale
        a = 0.33
        b = 10
        c = 0.33
        d = 1
        e = 0.33
        f = 0.1

        x = np.linspace(0,survival_cut*time_output_conv,survival_cut)
        y = triple_exp_func(x, a, b, c, d, e, f)
        #y = double_exp_func(x, a, b, c, d)
        yn = np.array(histo)

        # Perform the curve fitting from my histogram to mathematical function
        popt, pcov = curve_fit(triple_exp_func, x, yn)
        #popt, pcov = curve_fit(double_exp_func, x, yn)

        yy = triple_exp_func(x, *popt)
        #yy = double_exp_func(x, *popt)

        # For each of the fitting states, print 'er out.
        with open(prefix+"_statefit"+str(state_id),"w") as out:
            for survival_time, population in enumerate(yy):
                out.write(str(survival_time*time_output_conv)+" "
                          +str(population)+"\n")

        # Though we could use linregress to obtain R^2, correlation suffices
        fit_correlation = np.corrcoef(histo, yy)[0][1]
        fit_values.append((list(popt),list(np.sqrt(pcov.diagonal())),
                           fit_correlation))

    return fit_values

if __name__ == '__main__':
    parser = ArgumentParser(
    description='Produces distributions of rotamer state survival times')

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

    print "Writing rotamer survival states for plotting", args.dividers
    print write_rotamer_state_survival(data_f_dunk, data_f_states,
                                       traj_col=args.traj_col)
