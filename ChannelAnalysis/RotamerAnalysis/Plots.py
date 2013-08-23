#!/usr/bin/python
###############################################################################
#
# This script takes processed data from RotamerAnalysis functions and prepares
# awesome plots using matplotlib. Unfortunately, it's highly specific
# to the system being studied and often requires serious tweaking of ranges
# and scaling (especially for histograms).
#
# By Chris Ing, 2013 for Python 2.7
#
###############################################################################
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib.patches import Rectangle
from numpy import array, transpose
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

# This is for plotting rotamer-rotamer 2d histograms.
def plot_rotamer_2dhistograms(rotamer_2dhisto, state_dividers, label_stats,
                               plot_title="Rotamer-Rotamer Histograms",
                               prefix=None):

    # Initialize the figure and compute the number of rows that will
    # be needed in the figure to capture all the data.
    fig = plt.figure()

    rot_histo = rotamer_2dhisto[0]
    rot_xs = rotamer_2dhisto[1]
    rot_ys = rotamer_2dhisto[2]

    # For labelling purposes, generate the same combinations that
    # will be used to title the plots
    ion_name_def = ["R","G","B","P"]

    ax = fig.add_subplot(1,1,1)
    extent = [rot_xs[0], rot_xs[-1], rot_ys[0], rot_ys[-1]]

    cax = ax.imshow(rot_histo.transpose(), extent=extent, origin="lower",
                    interpolation='nearest', cmap=mpl.cm.CMRmap)
    ax.autoscale(False)

    # This prints a divider line based on the inputted state_divider
    # variables.
    for divider_val in state_dividers:
        ax.plot(range(int(rot_xs[0]),int(rot_xs[-1])),
                [divider_val]*(rot_xs[-1]-rot_xs[0]),
                linewidth=2.0,color='k')

    # This prints text on the 2D histogram that indictes the percentage
    # data in that dihedral space defined between the dividers
    num_regions = len(state_dividers)+1

    # This iterates over the state divider limits (top and bottom)
    # in order to position text within those boundaries.
    all_dividers = [rot_ys[0]]+state_dividers+[rot_ys[-1]]
    for range_id, range_vals in enumerate(window(all_dividers)):
        y_position = ((range_vals[0]+range_vals[1])/2)/rot_ys[-1]
        plt.text(0.90, y_position,
         r"{:.1f}%".format(100*label_stats[range_id]),
         #weight="bold",
         ha='center', va='center', transform=ax.transAxes,
         bbox=dict(facecolor='white', alpha=0.5))

    ax.yaxis.set_major_locator(MultipleLocator(60))
    ax.yaxis.set_minor_locator(MultipleLocator(30))
    ax.xaxis.set_major_locator(MultipleLocator(60))
    ax.xaxis.set_minor_locator(MultipleLocator(30))
    ax.set_axisbelow(True)
    ax.yaxis.grid(True,'minor')
    ax.yaxis.grid(True,'major', linewidth=0.5, linestyle='-')
    ax.xaxis.grid(True,'minor')
    ax.xaxis.grid(True,'major', linewidth=0.5, linestyle='-')

    ax.set_xlabel("Chi1 (degrees)")
    ax.set_ylabel("Chi2 (degrees)")

    #ax.set_title(plot_title)
    plt.subplots_adjust(hspace = 0.2, wspace = 0.02,
                        left = 0.1, bottom = 0.08,
                        right = 0.8, top = 0.90)
    plt.suptitle(plot_title)

    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.colorbar(cax, cax=cbar_ax, ticks=[0, 1, 2, 3, 4, 5, 6, 7])
    #cbar_ax.set_yticklabels(['free energy (kcal/mol)'], rotation=90)
    cbar_ax.set_ylabel('free energy (kcal/mol)', rotation=270)

    if prefix != None:
        plt.savefig(prefix+".pdf")
    else:
        plt.show()

    return True