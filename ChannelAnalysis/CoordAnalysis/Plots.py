#!/usr/bin/python
###############################################################################
#
# This script takes processed data from CoordAnalysis functions and prepares
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
from numpy import array
from itertools import combinations,product
from operator import itemgetter

# Plots ion binding modes along channel axis.
def plot_bindingmode_histograms(group_histo_1st,
                                group_histo_2nd,
                                group_histo_both,
                                selectivityf_map,
                                plot_title="Ion Binding Mode Histograms",
                                prefix=None):

    # Make a list of all binding modes so that the colors will be the same
    # in all of the plots.
    all_modes = sorted(set(group_histo_1st[0].keys() +
                           group_histo_2nd[0].keys() +
                           group_histo_both[0].keys()),reverse=True)

    # Labelling is very complicated for multiple subplots with different
    # lines on each one. All must be tracked.
    leg_series = []
    leg_lbl = []
    leg_num = []

    fig = plt.figure()

    # PLOT 1 - 1ST COORDINATION

    ax1 = fig.add_subplot(3,1,1)
    # 12 plot color scheme generated using http://colorbrewer2.org/
    # with extra ones added on the end...
    plot_colors = ["#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3",
                  "#FDB462","#B3DE69","#FCCDE5","#D9D9D9","#BC80BD",
                  "#CCEBC5","#FFED6F","#E31A1C","#6A3D9A","#1F78B4",
                  "#33A02C","#BF812D","#01665E","#C51B7D","#276419",
                  "#C51B7D","#BABABA","#8E0152","#FF7F00","#FFFF33"]

    num_sfvals = len(group_histo_1st[0].keys()[0])

    for mode_id, mode in enumerate(all_modes):
        if mode in group_histo_1st[0]:
            mode_histo = group_histo_1st[0][mode]
            mode_x = group_histo_1st[1][mode][1:]
            diff = (mode_x[1]-mode_x[0])/2.0
            mode_x_shift = (mode_x - diff)*-10

            # Here we treat the ++ (any coordinated ion) and 00
            # (no coordination) ions specifically.
            if mode == "+"*num_sfvals:
                mode_string = "Bound to SF"
                temp_plot, = ax1.plot(mode_x_shift, mode_histo, linewidth=2.5,
                                     color="black", zorder = 0,
                                     label="Bound to SF")
            elif mode == "0"*num_sfvals:
                mode_string = "Unbound to SF"
                temp_plot, = ax1.plot(mode_x_shift, mode_histo, 'k--', linewidth=2.5,
                                     color="black", zorder = 0,
                                     label="Unfound to SF")
            else:
                mode_string = " ".join([selectivityf_map[int(num)]
                                             for num, digit in enumerate(mode)
                                             if int(digit) > 0])
                #temp_plot, = ax1.plot(mode_x_shift, mode_histo,
                #                     linewidth=2.0,
                #                     color=plot_colors[mode_id],
                #                     label=mode_string)

                temp_plot = Rectangle((0.5, 0.5), 1, 1,
                                      color=plot_colors[mode_id], alpha=0.75)
                ax1.fill_between(mode_x_shift, 0, mode_histo,
                                             linewidth=0.5,
                                             facecolor=plot_colors[mode_id],
                                             alpha=0.75)

            if mode_string not in leg_lbl:
                leg_series.append(temp_plot)
                leg_lbl.append(mode_string)
                leg_num.append(mode)

    ax1.set_ylim(bottom=0)
    ax1.xaxis.set_major_locator(MultipleLocator(2))
    ax1.xaxis.set_minor_locator(MultipleLocator(1))
    #ax1.xaxis.set_major_locator(MultipleLocator(50))
    #ax1.xaxis.set_minor_locator(MultipleLocator(25))
    ax1.set_axisbelow(True)
    ax1.yaxis.grid(True,'minor')
    ax1.yaxis.grid(True,'major', linewidth=0.5, linestyle='-')
    ax1.xaxis.grid(True,'minor')
    ax1.xaxis.grid(True,'major', linewidth=0.5, linestyle='-')
    ax1.set_ylabel("1st")
    plt.setp(ax1.get_yticklabels(), visible=False)
    plt.setp(ax1.get_xticklabels(), visible=False)
    ax1.grid(True)

    # PLOT 2 - 2ND COORDINATION

    ax2 = fig.add_subplot(3,1,2)

    for mode_id, mode in enumerate(all_modes):
        if mode in group_histo_2nd[0]:
            mode_histo = group_histo_2nd[0][mode]
            mode_x = group_histo_2nd[1][mode][1:]
            diff = (mode_x[1]-mode_x[0])/2.0
            mode_x_shift = (mode_x - diff)*-10

            # Here we treat the ++ (any coordinated ion) and 00
            # (no coordination) ions specifically.
            if mode == "+"*num_sfvals:
                mode_string = "Bound to SF"
                temp_plot, = ax2.plot(mode_x_shift, mode_histo, linewidth=2.5,
                                     color="black", zorder = 0,
                                     label="Bound to SF")
            elif mode == "0"*num_sfvals:
                mode_string = "Unbound to SF"
                temp_plot, = ax2.plot(mode_x_shift, mode_histo, 'k--', linewidth=2.5,
                                     color="black", zorder = 0,
                                     label="Unfound to SF")
            else:
                mode_string = " ".join([selectivityf_map[int(num)]
                                             for num, digit in enumerate(mode)
                                             if int(digit) > 0])
                #temp_plot, = ax2.plot(mode_x_shift, mode_histo,
                #                     linewidth=2.0,
                #                     color=plot_colors[mode_id],
                #                     label=mode_string)

                temp_plot = Rectangle((0.5, 0.5), 1, 1,
                                      color=plot_colors[mode_id], alpha=0.75)
                ax2.fill_between(mode_x_shift, 0, mode_histo,
                                             linewidth=0.5,
                                             facecolor=plot_colors[mode_id],
                                             alpha=0.75)

            if mode_string not in leg_lbl:
                leg_series.append(temp_plot)
                leg_lbl.append(mode_string)
                leg_num.append(mode)

    ax2.set_ylim(bottom=0)
    ax2.xaxis.set_major_locator(MultipleLocator(2))
    ax2.xaxis.set_minor_locator(MultipleLocator(1))
    #ax2.xaxis.set_major_locator(MultipleLocator(50))
    #ax2.xaxis.set_minor_locator(MultipleLocator(25))
    ax2.set_axisbelow(True)
    ax2.yaxis.grid(True,'minor')
    ax2.yaxis.grid(True,'major', linewidth=0.5, linestyle='-')
    ax2.xaxis.grid(True,'minor')
    ax2.xaxis.grid(True,'major', linewidth=0.5, linestyle='-')
    ax2.set_ylabel("2nd")
    plt.setp(ax2.get_yticklabels(), visible=False)
    plt.setp(ax2.get_xticklabels(), visible=False)
    ax2.grid(True)

    # PLOT 3 - BOTH COORDINATION

    ax3 = fig.add_subplot(3,1,3)

    for mode_id, mode in enumerate(all_modes):
        if mode in group_histo_both[0]:
            mode_histo = group_histo_both[0][mode]
            mode_x = group_histo_both[1][mode][1:]
            diff = (mode_x[1]-mode_x[0])/2.0
            mode_x_shift = (mode_x - diff)*-10

            # Here we treat the ++ (any coordinated ion) and 00
            # (no coordination) ions specifically.
            if mode == "+"*num_sfvals:
                mode_string = "Bound to SF"
                temp_plot, = ax3.plot(mode_x_shift, mode_histo, linewidth=2.5,
                                     color="black", zorder = 0,
                                     label="Bound to SF")
            elif mode == "0"*num_sfvals:
                mode_string = "Unbound to SF"
                temp_plot, = ax3.plot(mode_x_shift, mode_histo, 'k--', linewidth=2.5,
                                     color="black", zorder = 0,
                                     label="Unfound to SF")
            else:
                mode_string = " ".join([selectivityf_map[int(num)]
                                             for num, digit in enumerate(mode)
                                             if int(digit) > 0])
                #temp_plot, = ax3.plot(mode_x_shift, mode_histo,
                #                     linewidth=2.0,
                #                     color=plot_colors[mode_id],
                #                     label=mode_string)

                temp_plot = Rectangle((0.5, 0.5), 1, 1,
                                      color=plot_colors[mode_id], alpha=0.75)
                ax3.fill_between(mode_x_shift, 0, mode_histo,
                                             linewidth=0.5,
                                             facecolor=plot_colors[mode_id],
                                             alpha=0.75)

            if mode_string not in leg_lbl:
                leg_series.append(temp_plot)
                leg_lbl.append(mode_string)
                leg_num.append(mode)

    ax3.set_ylim(bottom=0)
    ax3.xaxis.set_major_locator(MultipleLocator(2))
    ax3.xaxis.set_minor_locator(MultipleLocator(1))
    #ax3.xaxis.set_major_locator(MultipleLocator(50))
    #ax3.xaxis.set_minor_locator(MultipleLocator(25))
    ax3.set_axisbelow(True)
    ax3.yaxis.grid(True,'minor')
    ax3.yaxis.grid(True,'major', linewidth=0.5, linestyle='-')
    ax3.xaxis.grid(True,'minor')
    ax3.xaxis.grid(True,'major', linewidth=0.5, linestyle='-')
    ax3.set_ylabel("Both")
    plt.setp(ax3.get_yticklabels(), visible=False)
    ax3.grid(True)

    # Now general plot things

    all_modes_names = []
    for mode in all_modes:
        if ((mode in group_histo_1st[0]) or
            (mode in group_histo_2nd[0]) or
            (mode in group_histo_both[0])):
            if mode == "+"*num_sfvals:
                all_modes_names.append("Bound to SF")
            elif mode == "0"*num_sfvals:
                all_modes_names.append("Unbound to SF")
            else:
                all_modes_names.append(" ".join([selectivityf_map[int(num)]
                                          for num, digit in enumerate(mode)
                                          if int(digit) > 0]))

    # Sort the axis labels based on the numerical code sorted in reverse.
    # This is so the "Bound to SF" and "Unbound to SF" envelopes are last.
    sorted_legenddata = sorted(zip(leg_series,leg_lbl, leg_num),
                               key=itemgetter(2), reverse=True)
    sorted_series, sorted_lbl, _ = zip(*sorted_legenddata)
    print sorted_lbl

    leg = ax1.legend(tuple(sorted_series),tuple(sorted_lbl),
                     bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                     ncol=6, mode="expand", borderaxespad=0.)
    for t in leg.get_texts():
        t.set_fontsize(6)

    plt.subplots_adjust(hspace = 0.1, wspace = 0.02,
                        left = 0.1, bottom = 0.08,
                        right = 0.95, top = 0.75)
    plt.suptitle(plot_title)

    if prefix != None:
        plt.savefig(prefix+".pdf")
    else:
        plt.show()

    return True

# This is for plotting ion-ion 2d histograms. It is useful for visualizing
# correlations between ion positions for different occupancy states.
# Expects quite a specific data format of keys in the form "occ_id:ion1-ion2"
def plot_ionsplit_2dhistograms(ionsplit_2dhisto, ion_num_cutoff=3,
                               plot_title="Ion Split 2D Histograms",
                               prefix=None):

    # Initialize the figure and compute the number of rows that will
    # be needed in the figure to capture all the data.
    fig = plt.figure()

    num_rows = len(ionsplit_2dhisto[0].keys())
    num_cols = max([len(val) for key,val in ionsplit_2dhisto[0].iteritems()])
    print "Detected ", num_rows, " rows and ", num_cols, " cols"

    print ionsplit_2dhisto[0].keys(), [len(val) for key,val in ionsplit_2dhisto[0].iteritems()]

    for row_count, occ_id in enumerate(sorted(ionsplit_2dhisto[0].keys())):
        occ_histos = ionsplit_2dhisto[0][occ_id]
        occ_xs = ionsplit_2dhisto[1][occ_id]
        occ_ys = ionsplit_2dhisto[2][occ_id]

        # For labelling purposes, generate the same combinations that
        # will be used to title the plots
        ion_name_def = ["R","G","B","P"]
        occ_pairs = list(combinations(range(min(occ_id,ion_num_cutoff)),2))
        occ_pairs_name = [(ion_name_def[x],ion_name_def[y]) for x,y in occ_pairs]

        plot_position = 1+row_count*num_cols
        first_col = True
        for pair_id, (histo, x, y) in enumerate(zip(occ_histos, occ_xs, occ_ys)):

            extent = [-10*x[0], -10*x[-1], 10*y[-1], 10*y[0]]
            ax = fig.add_subplot(num_rows, num_cols, plot_position)
            histo_y_flipped = histo[::-1, :]

            cax = ax.imshow(histo_y_flipped, extent=extent,
                            interpolation='nearest', cmap=mpl.cm.CMRmap)
            ax.autoscale(False)
            ax.plot(range(-13,14),range(-13,14),linewidth=1.0,color='k')
            ax.yaxis.set_major_locator(MultipleLocator(4))
            ax.yaxis.set_minor_locator(MultipleLocator(2))
            ax.xaxis.set_major_locator(MultipleLocator(4))
            ax.xaxis.set_minor_locator(MultipleLocator(2))
            ax.set_axisbelow(True)
            ax.yaxis.grid(True,'minor')
            ax.yaxis.grid(True,'major', linewidth=0.5, linestyle='-')
            ax.xaxis.grid(True,'minor')
            ax.xaxis.grid(True,'major', linewidth=0.5, linestyle='-')
            ax.set_title(str(occ_pairs_name[pair_id]))
            plot_position += 1

            if row_count+1 != num_rows:
                plt.setp(ax.get_xticklabels(), visible=False)

            if not first_col:
                plt.setp(ax.get_yticklabels(), visible=False)
            else:
                ax.set_ylabel(str(occ_id)+" ions")

            first_col = False

    #ax.set_title(plot_title)
    plt.subplots_adjust(hspace = 0.2, wspace = 0.02,
                        left = 0.1, bottom = 0.08,
                        right = 0.8, top = 0.90)
    plt.suptitle(plot_title)
    # Add colorbar, make sure to specify tick locations to match desired ticklabels
    #plt.colorbar()
    #cbar = fig.colorbar(cax, ticks=[0, 2.5, 5])
    #fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.colorbar(cax, cax=cbar_ax, ticks=[0, 1, 2, 3, 4, 5, 6])
    #cbar_ax.set_yticklabels(['free energy (kcal/mol)'], rotation=90)
    cbar_ax.set_ylabel('free energy (kcal/mol)', rotation=270)


    if prefix != None:
        plt.savefig(prefix+".pdf")
    else:
        plt.show()

    return True

# Plots inner ion species ordering histograms.
def plot_species_histograms(species_histo, species_map, species_populations,
                            plot_title="Ion Species Ordering Histograms",
                            prefix=None):

    # Initialize the figure and compute the number of rows that will
    # be needed in the figure to capture all the data.
    fig = plt.figure()

    num_rows = len(species_histo[0].keys())

    for row_count, species_id in enumerate(sorted(species_histo[0].keys())):
        spc_histos =  species_histo[0][species_id]
        spc_xs = species_histo[1][species_id]

        # For labelling purposes, generate the same combinations that
        # will be used to title the plots
        color_list = ["R","G","B","#663399","C","0.1","0.2","0.3",
                      "0.4","0.5","0.6","0.7"]

        # First, extract the population of the particular occupancy
        # state, if it doesn't mean the cutoff of 0.5%< of the dataset
        # then give it a boot!
        spc_index = species_populations[1]["MEAN"].index(int(species_id))
        mean_val = species_populations[0]["MEAN"][spc_index]*100
        stderr_val = species_populations[0]["STDERR"][spc_index]*100

        ax = fig.add_subplot(num_rows, 1, row_count+1)
        for ion_id, (histo, edges) in enumerate(zip(spc_histos, spc_xs)):

            channel_spc_hist_x = array(edges[1:])
            diff = (channel_spc_hist_x[1]-channel_spc_hist_x[0])/2.0
            channel_spc_hist_x_shift = (channel_spc_hist_x - diff)*-10
            channel_spc_hist_y = array(histo)
            ax.plot(channel_spc_hist_x_shift, channel_spc_hist_y,
                    linewidth=2.5, color=color_list[ion_id])

        #ax.yaxis.set_major_locator(MultipleLocator(4))
        #ax.yaxis.set_minor_locator(MultipleLocator(2))
        ax.xaxis.set_major_locator(MultipleLocator(2))
        ax.xaxis.set_minor_locator(MultipleLocator(1))
        ax.set_axisbelow(True)
        ax.yaxis.grid(True,'minor')
        ax.yaxis.grid(True,'major', linewidth=0.5, linestyle='-')
        ax.xaxis.grid(True,'minor')
        ax.xaxis.grid(True,'major', linewidth=0.5, linestyle='-')
        ax.set_ylabel(species_map[int(species_id)])
        plt.setp(ax.get_yticklabels(), visible=False)
        plt.setp(ax.get_xticklabels(), visible=False)
        #l = ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        #for t in l.get_texts():
        #    t.set_fontsize(6)

        plt.text(0.15, 0.70,
                 r"Pop: {:.1f}% $\pm$ {:.1f}%".format(mean_val,stderr_val),
                 ha='center', va='center', transform=ax.transAxes)

    else:
        # If it's the last plot of the for loop, then give it special
        # treatment. Technically this could be outside the for clause...
        plt.setp(ax.get_xticklabels(), visible=True)
        ax.set_xlabel("Axial position (Ang)")

    plt.subplots_adjust(hspace = 0.1, wspace = 0.02,
                        left = 0.1, bottom = 0.08,
                        right = 0.80, top = 0.93)
    plt.suptitle(plot_title)

    if prefix != None:
        plt.savefig(prefix+".pdf")
    else:
        plt.show()

    return True

# Plots ion pair distance histograms.
def plot_iondist_histograms(iondist_histo, ion_num_cutoff=3,
                               plot_title="Ion Pair Distance Histograms",
                               prefix=None):

    # Initialize the figure and compute the number of rows that will
    # be needed in the figure to capture all the data.
    fig = plt.figure()

    num_rows = len(iondist_histo[0].keys())

    for row_count, occ_id in enumerate(sorted(iondist_histo[0].keys())):
        occ_histos = iondist_histo[0][occ_id]
        occ_xs = iondist_histo[1][occ_id]

        # For labelling purposes, generate the same combinations that
        # will be used to title the plots
        ion_name_def = ["R","G","B","P"]
        occ_pairs = list(combinations(range(min(occ_id,ion_num_cutoff)),2))
        occ_pairs_name = [(ion_name_def[x],ion_name_def[y]) for x,y in occ_pairs]

        ax = fig.add_subplot(num_rows, 1, row_count+1)
        for pair_id, (histo, edges) in enumerate(zip(occ_histos, occ_xs)):

            channel_occ_hist_x = array(edges[1:])
            diff = (channel_occ_hist_x[1]-channel_occ_hist_x[0])/2.0
            channel_occ_hist_x_shift = (channel_occ_hist_x - diff)*-10
            channel_occ_hist_y = array(histo)
            ax.plot(channel_occ_hist_x_shift, channel_occ_hist_y,
                    linewidth=2.5, label=str(occ_pairs_name[pair_id]))

        #ax.yaxis.set_major_locator(MultipleLocator(4))
        #ax.yaxis.set_minor_locator(MultipleLocator(2))
        ax.xaxis.set_major_locator(MultipleLocator(2))
        ax.xaxis.set_minor_locator(MultipleLocator(1))
        ax.set_axisbelow(True)
        ax.yaxis.grid(True,'minor')
        ax.yaxis.grid(True,'major', linewidth=0.5, linestyle='-')
        ax.xaxis.grid(True,'minor')
        ax.xaxis.grid(True,'major', linewidth=0.5, linestyle='-')
        ax.set_ylabel(str(occ_id)+" Ion Pair Dens.")
        plt.setp(ax.get_yticklabels(), visible=False)
        plt.setp(ax.get_xticklabels(), visible=False)
        l = ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        for t in l.get_texts():
            t.set_fontsize(6)

    else:
        # If it's the last plot of the for loop, then give it special
        # treatment. Technically this could be outside the for clause...
        plt.setp(ax.get_xticklabels(), visible=True)
        ax.set_xlabel("Axial position (Ang)")

    plt.subplots_adjust(hspace = 0.1, wspace = 0.02,
                        left = 0.1, bottom = 0.08,
                        right = 0.80, top = 0.93)
    plt.suptitle(plot_title)

    if prefix != None:
        plt.savefig(prefix+".pdf")
    else:
        plt.show()

    return True

# This is for plotting global ion distributions for different
# ion occupancy states along with statistics about how often
# those states occured in the entire dataset. This is a figure of
# the global data and shows averages across all trajectories.
def plot_ionsplit_histograms(ionsplit_histo, occ_populations,
                             allion_histo, allatom_histos,
                             occ_pop_cutoff = 0.5,
                             plot_title="Ion Split Histograms",
                             prefix=None):

    # Color cycle should be initialized very early.
    mpl.rcParams['axes.color_cycle'] = ['r', 'g', 'b', 'm', 'c', 'y', 'k']

    # Initialize the figure and compute the number of rows that will
    # be needed in the figure to capture all the data.
    fig = plt.figure()

    # Prepare a list of binary values and use it to compute how many
    # are above the occupancy population cutoff.
    num_plots = sum([x*100 >= occ_pop_cutoff
                    for x in occ_populations[0]["MEAN"]])+1

    # Initialize the top-most histogram
    ax = fig.add_subplot(num_plots,1,1)

    # Define a new color scheme for this purpose only.
    oxy_colors = ["#DCFFAC","#E8AB4E","#FF49A4","#4E63E8","#56FFA6",
                  "#248077","#74AD8D","#C82754","#F7BB21","#F9E2B7"]

    # Plot the channel oxygen distributions, however many there are.
    for color_id, histo in enumerate(allatom_histos):
        #for color_id, histo in enumerate(allatom_histos[::-1]):
        channel_oxy_hist_x = (array(histo[1]["ALL"][1:]))
        diff = (channel_oxy_hist_x[1]-channel_oxy_hist_x[0])/2.0
        channel_oxy_hist_x_shift = (channel_oxy_hist_x - diff)*-1
        channel_oxy_hist_y = (array(histo[0]["ALL"]))
        ax.fill_between(channel_oxy_hist_x_shift, 0, channel_oxy_hist_y,
                linewidth=0.5, facecolor=oxy_colors[color_id], alpha=0.8)

    channel_occ_hist_x = (array(allion_histo[1]["ALL"][1:]))
    diff = (channel_occ_hist_x[1]-channel_occ_hist_x[0])/2.0
    channel_occ_hist_x_shift = (channel_occ_hist_x - diff)*-1

    # Here we multiply the height by 2 since typically it isn't as high
    # as the oxygen distribution plotted above.
    channel_occ_hist_y = (array(allion_histo[0]["ALL"]))*2
    ax.plot(channel_occ_hist_x_shift, channel_occ_hist_y,
            linewidth=2.5, color="black")

    ax.set_ylim(bottom=0)
    ax.grid(True)
    ax.set_ylabel("ALL Density")
    plt.setp(ax.get_xticklabels(), visible=False)
    plt.setp(ax.get_yticklabels(), visible=False)
    ax.xaxis.set_major_locator(MultipleLocator(0.2))
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))

    # Loop over the datafile which is a tuple of two dictionaries
    # (one with the z vals and the other with the histogram height)
    # that are both key'd by the occupancy number and valued by
    # a list of length equal to that occupancy number.
    occ_count = 1
    for occ_id in ionsplit_histo[0].keys():

        # First, extract the population of the particular occupancy
        # state, if it doesn't mean the cutoff of 0.5%< of the dataset
        # then give it a boot!
        occ_index = occ_populations[1]["MEAN"].index(occ_id)
        mean_val = occ_populations[0]["MEAN"][occ_index]*100
        stderr_val = occ_populations[0]["STDERR"][occ_index]*100

        if mean_val >= occ_pop_cutoff:

            # Since it's only 1 column of a known number of plots,
            # add it to the existing figure increasingly lower.
            ax = fig.add_subplot(num_plots,1,occ_count+1)

            # Increase the counter of valid occupancy states
            occ_count += 1

            # Now plot multiple lines on each plot depending on the
            # length of the occupancy array.
            for ion_num in range(len(ionsplit_histo[1][occ_id])):
                channel_occ_hist_x = (array(ionsplit_histo[1][occ_id][ion_num][1:]))
                diff = (channel_occ_hist_x[1]-channel_occ_hist_x[0])/2.0
                channel_occ_hist_x_shift = (channel_occ_hist_x - diff)*-10
                channel_occ_hist_y = (array(ionsplit_histo[0][occ_id][ion_num]))
                ax.plot(channel_occ_hist_x_shift, channel_occ_hist_y, linewidth=2.5)

            ax.grid(True)
            ax.set_ylim(bottom=0)
            ax.set_ylabel(str(occ_id)+" Ion Density")
            ax.set_xlabel("Axial position (Ang)")
            plt.setp(ax.get_xticklabels(), visible=False)
            plt.setp(ax.get_yticklabels(), visible=False)
            ax.xaxis.set_major_locator(MultipleLocator(2))
            ax.xaxis.set_minor_locator(MultipleLocator(1))

            plt.text(0.15, 0.9,
                     r"Pop: {:.1f}% $\pm$ {:.1f}%".format(mean_val,stderr_val),
                     ha='center', va='center', transform=ax.transAxes)
    else:
        # If it's the last plot of the for loop, then give it special
        # treatment. Technically this could be outside the for clause...
        plt.setp(ax.get_xticklabels(), visible=True)
        ax.set_xlabel("Axial position (nm)")

    plt.subplots_adjust(hspace = 0.1, wspace = 0.02,
                        left = 0.1, bottom = 0.08,
                        right = 0.95, top = 0.93)
    plt.suptitle(plot_title)

    if prefix != None:
        plt.savefig(prefix+".pdf")
    else:
        plt.show()

    return True

# This plots a stacked timeseries of trajectory properties in the main column
# and histograms of this data in a second column. Sorry for the god-awful
# number of arguments...
def plot_oxy_timeseries(alloxygen_ts, sf_ts,
                        time_conv=0.02,
                        prefix=None,
                        plot_title="Oxygen Timeseries",
                        max_coord=4,
                        max_length=500,
                        data_skip=10,
                        hist_scale=0.4):

    # This iterates over all the trajectory id's that you computed data for.
    # ion_timeseries[0] is the time values array, but any index would suffice.
    for traj_id in alloxygen_ts[0][0].keys():

        # Color cycle should be initialized very early.
        mpl.rcParams['axes.color_cycle'] = ['r', 'g', 'b', 'm', 'c', 'y', 'k']

        # Initialize the figure and compute the number of rows that will
        # be needed in the figure to capture all the data.
        fig = plt.figure()
        gs = gridspec.GridSpec(2, 1, height_ratios=[4,1])
        ax1 = plt.subplot(gs[0])

        # Prepare a list of binary values and use it to compute how many
        # are above the occupancy population cutoff.
        num_plots = len(alloxygen_ts)

        for row_id, oxygens in enumerate(alloxygen_ts):

            # Since it's only 1 column of a known number of plots,
            # add it to the existing figure increasingly lower.


            # Plot 1 - Ion Timeseries
            colors = ['r', 'g', 'b', 'm', 'c', 'y', 'k']
            for color_id, atom_id in enumerate(sorted(oxygens[0][traj_id].keys())):
                atom_ts_x = oxygens[1][traj_id][atom_id][::int(data_skip)]
                atom_ts_y = oxygens[0][traj_id][atom_id][::int(data_skip)]
                atom_ts_x_scale = [time_conv*pos for pos in atom_ts_x]
                atom_ts_y_flip = [-10*pos for pos in atom_ts_y]
                ax1.scatter(atom_ts_x_scale, atom_ts_y_flip, s=0.5,
                           color=colors[color_id])

        ax1.set_ylabel("Axial position (nm)")
        #ax1.set_ylim([-8, 12])
        ax1.set_ylim([-5, 7])
        #ax1.set_yticks([x/1.0 for x in range(-8,13,2)])
        ax1.set_yticks([x/1.0 for x in range(-5,8,1)])
        ax1.set_xlim([-5,max_length+5])
        ax1.set_xticks(range(0, int(50*round(max_length/50)), 50))
        ax1.yaxis.set_major_locator(MultipleLocator(2))
        ax1.yaxis.set_minor_locator(MultipleLocator(1))
        ax1.xaxis.set_major_locator(MultipleLocator(50))
        ax1.xaxis.set_minor_locator(MultipleLocator(25))
        ax1.set_axisbelow(True)
        ax1.yaxis.grid(True,'minor')
        ax1.yaxis.grid(True,'major', linewidth=0.5, linestyle='-')
        ax1.xaxis.grid(True,'minor')
        ax1.xaxis.grid(True,'major', linewidth=0.5, linestyle='-')
        plt.setp(ax1.get_xticklabels(), visible=False)

        ax2 = plt.subplot(gs[1], sharex=ax1)
        for color_id, atom_id in enumerate(sorted(sf_ts[0][traj_id].keys())):
            atom_ts_x = sf_ts[1][traj_id][atom_id][::int(data_skip)]
            atom_ts_y = sf_ts[0][traj_id][atom_id][::int(data_skip)]
            atom_ts_x_scale = [time_conv*pos for pos in atom_ts_x]
            # Awkwardly, the units for these files are incorrect and should
            # be fixed. They are in Ang already.
            atom_ts_y_flip = [-1*pos for pos in atom_ts_y]
            ax2.scatter(atom_ts_x_scale, atom_ts_y_flip, s=0.5,
                       color=colors[color_id])

        ax2.set_ylabel("C.O.M. Dev. (Ang)")
        ax2.set_xlabel("Time (ns)")
        ax2.set_ylim([-5, 5])
        ax1.set_xlim([-5,max_length+5])
        ax1.set_xticks(range(0, int(50*round(max_length/50)), 50))
        ax2.set_yticks([x/1.0 for x in range(-5,6,2)])
        ax2.yaxis.set_major_locator(MultipleLocator(1.0))
        ax2.yaxis.set_minor_locator(MultipleLocator(0.5))
        ax2.yaxis.grid(True,'minor')
        ax2.yaxis.grid(True,'major', linewidth=0.5, linestyle='-')
        ax2.xaxis.grid(True,'minor')
        ax2.xaxis.grid(True,'major', linewidth=0.5, linestyle='-')


        plt.subplots_adjust(hspace = 0.1, wspace = 0.02,
                            left = 0.1, bottom = 0.08,
                            right = 0.95, top = 0.93)
        plt.suptitle(plot_title+" - Traj "+str(traj_id))

        if prefix != None:
            plt.savefig(prefix+"_"+str(traj_id)+".pdf")
        else:
            plt.show()

        plt.close()

    return True

# Just a little helper method to attach labels
def autolabel(rects, ax):
    # attach some text labels
    for rect in rects:
        height = rect.get_height()
        ax.text(rect.get_x()+rect.get_width()/2., 1.05*height, '%d'% int(100*height),
                ha='center', va='bottom')

# This prints each binding mode and the mean occupancy of each of them
# across the dataset. Future support will be for multiple atom species.
def plot_mode_statistics(mode_counts, sf_col, selectivityf_code,
                         prefix=None, plot_title="Binding Mode Statistics"):

    num_plots = len(mode_counts)

    # Compute all the binding modes and turn them into strings
    all_binding_modes = product("01",repeat=len(sf_col))
    all_binding_strs = ["".join(mode) for mode in all_binding_modes]
    all_binding_chrstrs = []

    colors = ['r', 'g', 'b', 'm', 'c', 'y', 'k']

    for binding_str in all_binding_strs:
        temp_str = []
        temp_bool = True
        for coord_index, coord_res in enumerate(binding_str):
            if int(coord_res) > 0:
                temp_bool=False
                temp_str.append(selectivityf_code[coord_index])
        if temp_bool:
            temp_str.append("N/A")
        all_binding_chrstrs.append("".join(temp_str))

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)

    width = 0.35
    for ion_num, ion_species in enumerate(mode_counts):
        data_mean_y = ion_species[0]["MEAN"]
        data_mean_x = ion_species[1]["MEAN"]
        data_stderr = ion_species[0]["STDERR"]

        rects = ax.bar(array(data_mean_x)+width*ion_num, data_mean_y, width,
                        color=colors[ion_num], yerr=data_stderr)

        autolabel(rects, ax)


    ax.set_xticklabels( tuple(all_binding_chrstrs) )

    plt.subplots_adjust(hspace = 0.1, wspace = 0.02,
                        left = 0.1, bottom = 0.08,
                        right = 0.95, top = 0.93)
    plt.suptitle(plot_title+" - Mean Vals")

    if prefix != None:
        plt.savefig(prefix+".pdf")
    else:
        plt.show()

    plt.close()


# This plots a stacked timeseries of dihedral values over time. Useful
# for correlating with the oxygen timeseries above.
def plot_dihedral_timeseries(dihe_ts,
                             time_conv=0.02,
                             prefix=None,
                             plot_title="Dihedral Timeseries",
                             max_length=500,
                             data_skip=10):

    colors = ['r', 'g', 'b', 'm', 'c', 'y', 'k']

    # This iterates over all the trajectory id's that you computed data for.
    # ion_timeseries[0] is the time values array, but any index would suffice.
    for traj_id in dihe_ts[0].keys():

        # Initialize the figure and compute the number of rows that will
        # be needed in the figure to capture all the data.
        fig = plt.figure()

        # This should be 4 in the general case, not exactly sure though...
        #num_plots = len(dihe_ts[0][0][0])
        num_plots = 4

        for color_id, res_id in enumerate(sorted(dihe_ts[0][traj_id].keys())):
            ax = fig.add_subplot(num_plots,1,res_id+1)
            atom_ts_x = dihe_ts[1][traj_id][res_id][::int(data_skip)]
            atom_ts_x_scale = array([time_conv*pos for pos in atom_ts_x])
            atom_ts_y = dihe_ts[0][traj_id][res_id][::int(data_skip)]

            ax.scatter(atom_ts_x_scale, atom_ts_y, s=0.5, color=colors[color_id])
            ax.set_ylim([0, 360])
            ax.set_xlim([-5,max_length+5])
            ax.set_xticks(range(0, int(50*round(max_length/50)), 50))
            ax.yaxis.set_major_locator(MultipleLocator(60.0))
            ax.yaxis.set_minor_locator(MultipleLocator(30.0))
            ax.yaxis.grid(True,'minor')
            ax.yaxis.grid(True,'major', linewidth=0.5, linestyle='-')
            ax.xaxis.grid(True,'minor')
            ax.xaxis.grid(True,'major', linewidth=0.5, linestyle='-')
            ax.set_ylabel("Chain "+str(color_id+1)+ " Angle")
            plt.setp(ax.get_xticklabels(), visible=False)

        plt.setp(ax.get_xticklabels(), visible=True)
        ax.set_xlabel("Time (ns)")


        plt.subplots_adjust(hspace = 0.1, wspace = 0.02,
                            left = 0.1, bottom = 0.08,
                            right = 0.95, top = 0.93)
        plt.suptitle(plot_title+" - Traj "+str(traj_id))

        if prefix != None:
            plt.savefig(prefix+"_"+str(traj_id)+".pdf")
        else:
            plt.show()

        plt.close()

    return True

# This plots a stacked timeseries of trajectory properties in the main column
# and histograms of this data in a second column. Sorry for the god-awful
# number of arguments...
def plot_sf_w_2res_timeseries(channel_occ, channel_counts,
                              dunking_ts, dunking_counts,
                              ion_ts_1st, ion_histo_1st,
                              ion_ts_2nd, ion_histo_2nd,
                              ion_ts_both, ion_histo_both,
                              ion_ts_e_1st, ion_histo_e_1st,
                              ion_ts_l_1st, ion_histo_l_1st,
                              ion_ts_e_2nd, ion_histo_e_2nd,
                              ion_ts_l_2nd, ion_histo_l_2nd,
                              ion_ts_e_both, ion_histo_e_both,
                              ion_ts_l_both, ion_histo_l_both,
                              ion_ts_r, ion_histo_r,
                              time_conv=0.02,
                              prefix=None,
                              plot_title="Stacked Timeseries",
                              max_coord=4,
                              max_length=500,
                              data_skip=10,
                              hist_scale=0.4,
                              ion_types_ts=None,
                              coord_cols_map=["E177","L176"]):

    # This iterates over all the trajectory id's that you computed data for.
    # ion_timeseries[0] is the time values array, but any index would suffice.
    for traj_id in ion_ts_1st[0].keys():

        fig = plt.figure()
        # The grid is a 2 column figure, with 2 + 6 rows.
        gs = gridspec.GridSpec(9, 2,
                               width_ratios=[4,1],
                               height_ratios=[4,1,1,1,1,1,1,1,2])
        mpl.rcParams['axes.color_cycle'] = ['r', 'g', 'b', 'm', 'c', 'y', 'k']

        # Here we initialize all our subplots on the grid defined above.
        ax1 = plt.subplot(gs[0])
        ax2 = plt.subplot(gs[1], sharey = ax1)
        ax3 = plt.subplot(gs[2], sharex = ax1)
        ax4 = plt.subplot(gs[3], sharey = ax3, sharex = ax2)
        ax5 = plt.subplot(gs[4], sharex = ax1)
        ax6 = plt.subplot(gs[5], sharey = ax5, sharex = ax2)
        ax7 = plt.subplot(gs[6], sharex = ax1)
        ax8 = plt.subplot(gs[7], sharey = ax5, sharex = ax2)
        ax9 = plt.subplot(gs[8], sharex = ax1)
        ax10 = plt.subplot(gs[9], sharey = ax5, sharex = ax2)
        ax11 = plt.subplot(gs[10], sharex = ax1)
        ax12 = plt.subplot(gs[11], sharey = ax5, sharex = ax2)
        ax13 = plt.subplot(gs[12], sharex = ax1)
        ax14 = plt.subplot(gs[13], sharey = ax5, sharex = ax2)
        ax15 = plt.subplot(gs[14], sharex = ax1)
        ax16 = plt.subplot(gs[15], sharey = ax5, sharex = ax2)
        ax17 = plt.subplot(gs[16], sharex = ax1)
        ax18 = plt.subplot(gs[17], sharey = ax17, sharex = ax2)

        # Many axis are removed.
        plt.setp(ax1.get_xticklabels(), visible=False)
        plt.setp(ax2.get_xticklabels(), visible=False)
        plt.setp(ax2.get_yticklabels(), visible=False)
        plt.setp(ax3.get_xticklabels(), visible=False)
        plt.setp(ax4.get_yticklabels(), visible=False)
        plt.setp(ax4.get_xticklabels(), visible=False)
        plt.setp(ax5.get_xticklabels(), visible=False)
        plt.setp(ax6.get_yticklabels(), visible=False)
        plt.setp(ax6.get_xticklabels(), visible=False)
        plt.setp(ax7.get_xticklabels(), visible=False)
        plt.setp(ax8.get_yticklabels(), visible=False)
        plt.setp(ax8.get_xticklabels(), visible=False)
        plt.setp(ax9.get_xticklabels(), visible=False)
        plt.setp(ax10.get_yticklabels(), visible=False)
        plt.setp(ax10.get_xticklabels(), visible=False)
        plt.setp(ax11.get_xticklabels(), visible=False)
        plt.setp(ax12.get_yticklabels(), visible=False)
        plt.setp(ax12.get_xticklabels(), visible=False)
        plt.setp(ax13.get_xticklabels(), visible=False)
        plt.setp(ax14.get_yticklabels(), visible=False)
        plt.setp(ax14.get_xticklabels(), visible=False)
        plt.setp(ax15.get_xticklabels(), visible=False)
        plt.setp(ax16.get_yticklabels(), visible=False)
        plt.setp(ax16.get_xticklabels(), visible=False)
        plt.setp(ax18.get_xticklabels(), visible=False)
        plt.setp(ax18.get_yticklabels(), visible=False)

        # Plot 1 - Ion Timeseries
        colors = ['r', 'g', 'b', 'm', 'c', 'y', 'k']
        for color_id, ion_id in enumerate(sorted(ion_ts_1st[0][traj_id].keys(),
                                                 reverse=False)):
            ion_ts_x = ion_ts_1st[1][traj_id][ion_id][::int(data_skip/2)]
            ion_ts_y = ion_ts_1st[0][traj_id][ion_id][::int(data_skip/2)]
            ion_ts_x_scale = array([time_conv*pos for pos in ion_ts_x])
            ion_ts_y_flip = array([-1*pos for pos in ion_ts_y])


            # If ion types are passed into this array then we make a mask
            # of all the ions with type "1" as opposed to "0", and then
            # use that to extract all the points to replot them with a different
            # color.
            if ion_types_ts != None:
                mask_of_old_ions = array(ion_types_ts[0][traj_id][ion_id][::int(data_skip/2)])==0
                mask_of_new_ions = array(ion_types_ts[0][traj_id][ion_id][::int(data_skip/2)])==1
                ax1.scatter(ion_ts_x_scale[mask_of_old_ions],
                            ion_ts_y_flip[mask_of_old_ions],
                            s=1, color=colors[color_id])
                ax1.scatter(ion_ts_x_scale[mask_of_new_ions],
                            ion_ts_y_flip[mask_of_new_ions], marker='^',
                            s=1, color=colors[color_id])
            else:
                ax1.scatter(ion_ts_x_scale, ion_ts_y_flip, s=1, color=colors[color_id])

        # Plot 2 - Ion Position Histogram
        for vals, ion_z in zip(ion_histo_1st[0][traj_id], ion_histo_1st[1][traj_id]):
            mpl.rcParams['axes.color_cycle'] = ['r', 'g', 'b', 'm', 'c', 'y', 'k']
            diff = (ion_z[1]-ion_z[0])/2.0
            ion_z_shifted = [-1*(ion-diff) for ion in ion_z[1:]]
            vals_scaled = [val/10.0 for val in vals]
            ax2.plot(vals_scaled, ion_z_shifted)

        # Plot 3 - Dunking Counts
        dunking_ts_x = dunking_ts[1][traj_id][::10]
        dunking_ts_x_scale = [time_conv*pos for pos in dunking_ts_x]
        dunking_ts_y = dunking_ts[0][traj_id][::10]
        ax3.plot(dunking_ts_x_scale, dunking_ts_y)

        # Plot 4 - Channel Occupancy Histogram
        channel_occ_hist_x = dunking_counts[1][traj_id]
        channel_occ_hist_y = dunking_counts[0][traj_id]
        ax4.barh(channel_occ_hist_x, channel_occ_hist_y, align='center')

        # Plot 5 - E Coordination Timeseries Per Atom
        for ion_id in sorted(ion_ts_1st[0][traj_id].keys(), reverse=False):
            mpl.rcParams['axes.color_cycle'] = ['r', 'g', 'b', 'm', 'c', 'y', 'k']
            ion_ts_x = ion_ts_e_1st[1][traj_id][ion_id][::data_skip]
            ion_ts_y = ion_ts_e_1st[0][traj_id][ion_id][::data_skip]
            ion_ts_x_scale = [time_conv*pos for pos in ion_ts_x]
            ion_ts_y_flip = [1*pos for pos in ion_ts_y]
            ax5.plot(ion_ts_x_scale, ion_ts_y_flip)

        # Plot 6 - Coordination Histogram
        #mpl.rcParams['axes.color_cycle'] = ['r', 'g', 'b', 'm', 'c', 'y', 'k']
        left = 0
        for color_id, hist_data in enumerate(zip(ion_histo_e_1st[0][traj_id],
                                                      ion_histo_e_1st[1][traj_id])):
            vals_scaled = array([val*hist_scale for val in hist_data[0]])
            labels_scaled = array(hist_data[1][1:])-1
            ax6.barh(labels_scaled, vals_scaled, align='center',
                     color=colors[color_id], left=left)
            left += vals_scaled

        # Plot 7 - L Coordination Timeseries Per Atom
        for ion_id in sorted(ion_ts_1st[0][traj_id].keys(), reverse=False):
            mpl.rcParams['axes.color_cycle'] = ['r', 'g', 'b', 'm', 'c', 'y', 'k']
            ion_ts_x = ion_ts_l_1st[1][traj_id][ion_id][::data_skip]
            ion_ts_y = ion_ts_l_1st[0][traj_id][ion_id][::data_skip]
            ion_ts_x_scale = [time_conv*pos for pos in ion_ts_x]
            ion_ts_y_flip = [1*pos for pos in ion_ts_y]
            ax7.plot(ion_ts_x_scale, ion_ts_y_flip)

        # Plot 8 - Coordination Histogram
        #mpl.rcParams['axes.color_cycle'] = ['r', 'g', 'b', 'm', 'c', 'y', 'k']
        left = 0
        for color_id, hist_data in enumerate(zip(ion_histo_l_1st[0][traj_id],
                                                      ion_histo_l_1st[1][traj_id])):
            vals_scaled = array([val*hist_scale for val in hist_data[0]])
            labels_scaled = array(hist_data[1][1:])-1
            ax8.barh(labels_scaled, vals_scaled, align='center',
                     color=colors[color_id], left=left)
            left += vals_scaled

        # Plot 9 - E Coordination Timeseries Per Atom
        for ion_id in sorted(ion_ts_2nd[0][traj_id].keys(), reverse=False):
            mpl.rcParams['axes.color_cycle'] = ['r', 'g', 'b', 'm', 'c', 'y', 'k']
            ion_ts_x = ion_ts_e_2nd[1][traj_id][ion_id][::data_skip]
            ion_ts_y = ion_ts_e_2nd[0][traj_id][ion_id][::data_skip]
            ion_ts_x_scale = [time_conv*pos for pos in ion_ts_x]
            ion_ts_y_flip = [1*pos for pos in ion_ts_y]
            ax9.plot(ion_ts_x_scale, ion_ts_y_flip)

        # Plot 6 - Coordination Histogram
        #mpl.rcParams['axes.color_cycle'] = ['r', 'g', 'b', 'm', 'c', 'y', 'k']
        left = 0
        for color_id, hist_data in enumerate(zip(ion_histo_e_2nd[0][traj_id],
                                                      ion_histo_e_2nd[1][traj_id])):
            vals_scaled = array([val*hist_scale for val in hist_data[0]])
            labels_scaled = array(hist_data[1][1:])-1
            ax10.barh(labels_scaled, vals_scaled, align='center',
                     color=colors[color_id], left=left)
            left += vals_scaled

        # Plot 11 - L Coordination Timeseries Per Atom
        for ion_id in sorted(ion_ts_2nd[0][traj_id].keys(), reverse=False):
            mpl.rcParams['axes.color_cycle'] = ['r', 'g', 'b', 'm', 'c', 'y', 'k']
            ion_ts_x = ion_ts_l_2nd[1][traj_id][ion_id][::data_skip]
            ion_ts_y = ion_ts_l_2nd[0][traj_id][ion_id][::data_skip]
            ion_ts_x_scale = [time_conv*pos for pos in ion_ts_x]
            ion_ts_y_flip = [1*pos for pos in ion_ts_y]
            ax11.plot(ion_ts_x_scale, ion_ts_y_flip)

        # Plot 6 - Coordination Histogram
        #mpl.rcParams['axes.color_cycle'] = ['r', 'g', 'b', 'm', 'c', 'y', 'k']
        left = 0
        for color_id, hist_data in enumerate(zip(ion_histo_l_2nd[0][traj_id],
                                                      ion_histo_l_2nd[1][traj_id])):
            vals_scaled = array([val*hist_scale for val in hist_data[0]])
            labels_scaled = array(hist_data[1][1:])-1
            ax12.barh(labels_scaled, vals_scaled, align='center',
                     color=colors[color_id], left=left)
            left += vals_scaled

        # Plot 13 - E Coordination Timeseries Per Atom
        for ion_id in sorted(ion_ts_both[0][traj_id].keys(), reverse=False):
            mpl.rcParams['axes.color_cycle'] = ['r', 'g', 'b', 'm', 'c', 'y', 'k']
            ion_ts_x = ion_ts_e_both[1][traj_id][ion_id][::data_skip]
            ion_ts_y = ion_ts_e_both[0][traj_id][ion_id][::data_skip]
            ion_ts_x_scale = [time_conv*pos for pos in ion_ts_x]
            ion_ts_y_flip = [1*pos for pos in ion_ts_y]
            ax13.plot(ion_ts_x_scale, ion_ts_y_flip)

        # Plot 14 - Coordination Histogram
        #mpl.rcParams['axes.color_cycle'] = ['r', 'g', 'b', 'm', 'c', 'y', 'k']
        left = 0
        for color_id, hist_data in enumerate(zip(ion_histo_e_both[0][traj_id],
                                                      ion_histo_e_both[1][traj_id])):
            vals_scaled = array([val*hist_scale for val in hist_data[0]])
            labels_scaled = array(hist_data[1][1:])-1
            ax14.barh(labels_scaled, vals_scaled, align='center',
                     color=colors[color_id], left=left)
            left += vals_scaled

        # Plot 15 - L Coordination Timeseries Per Atom
        for ion_id in sorted(ion_ts_both[0][traj_id].keys(), reverse=False):
            mpl.rcParams['axes.color_cycle'] = ['r', 'g', 'b', 'm', 'c', 'y', 'k']
            ion_ts_x = ion_ts_l_both[1][traj_id][ion_id][::data_skip]
            ion_ts_y = ion_ts_l_both[0][traj_id][ion_id][::data_skip]
            ion_ts_x_scale = [time_conv*pos for pos in ion_ts_x]
            ion_ts_y_flip = [1*pos for pos in ion_ts_y]
            ax15.plot(ion_ts_x_scale, ion_ts_y_flip)

        # Plot 16 - Coordination Histogram
        #mpl.rcParams['axes.color_cycle'] = ['r', 'g', 'b', 'm', 'c', 'y', 'k']
        left = 0
        for color_id, hist_data in enumerate(zip(ion_histo_l_both[0][traj_id],
                                                      ion_histo_l_both[1][traj_id])):
            vals_scaled = array([val*hist_scale for val in hist_data[0]])
            labels_scaled = array(hist_data[1][1:])-1
            ax16.barh(labels_scaled, vals_scaled, align='center',
                     color=colors[color_id], left=left)
            left += vals_scaled

        # Plot 17 - Ion R Timeseries
        colors = ['r', 'g', 'b', 'm', 'c', 'y', 'k']
        for color_id, ion_id in enumerate(sorted(ion_ts_r[0][traj_id].keys(),
                                                 reverse=False)):
            ion_ts_x = ion_ts_r[1][traj_id][ion_id][::data_skip]
            ion_ts_y = ion_ts_r[0][traj_id][ion_id][::data_skip]
            ion_ts_x_scale = [time_conv*pos for pos in ion_ts_x]
            ion_ts_y_flip = [1*pos for pos in ion_ts_y]
            ax17.scatter(ion_ts_x_scale, ion_ts_y_flip, s=1, color=colors[color_id])

        # Plot 18 - Ion R Position Histogram
        for vals, ion_z in zip(ion_histo_r[0][traj_id], ion_histo_r[1][traj_id]):
            mpl.rcParams['axes.color_cycle'] = ['r', 'g', 'b', 'm', 'c', 'y', 'k']
            diff = (ion_z[1]-ion_z[0])/2.0
            ion_z_shifted = [1*(ion-diff) for ion in ion_z[1:]]
            vals_scaled = [val/20.0 for val in vals]
            ax18.plot(vals_scaled, ion_z_shifted)

        ax1.set_ylabel("Axial position (nm)")
        ax1.set_ylim([-1.1, 1.0])
        ax1.set_yticks([x/10.0 for x in range(-10,10,2)])
        ax1.set_xlim([-5,max_length+5])
        ax1.set_xticks(range(0, int(50*round(max_length/50)), 50))
        ax1.yaxis.set_major_locator(MultipleLocator(0.2))
        ax1.yaxis.set_minor_locator(MultipleLocator(0.1))
        ax1.xaxis.set_major_locator(MultipleLocator(50))
        ax1.xaxis.set_minor_locator(MultipleLocator(25))
        ax1.set_axisbelow(True)
        ax1.yaxis.grid(True,'minor')
        ax1.yaxis.grid(True,'major', linewidth=0.5, linestyle='-')
        ax1.xaxis.grid(True,'minor')
        ax1.xaxis.grid(True,'major', linewidth=0.5, linestyle='-')

        ax2.grid(True)
        ax2.set_xlim([0,1.0])
        ax2.set_xticks([x/10.0 for x in range(0,11,2)])

        ax3.set_ylabel("Dunk")
        ax3.set_yticks(range(0,max_coord+1,2))
        ax3.set_ylim([-0.5,max_coord+0.5])
        ax3.yaxis.set_major_locator(MultipleLocator(2))
        ax3.yaxis.set_minor_locator(MultipleLocator(1))
        ax3.xaxis.set_major_locator(MultipleLocator(50))
        ax3.xaxis.set_minor_locator(MultipleLocator(25))
        ax3.set_axisbelow(True)
        ax3.yaxis.grid(True,'minor')
        ax3.yaxis.grid(True,'major', linewidth=0.5, linestyle='-')
        ax3.xaxis.grid(True,'minor')
        ax3.xaxis.grid(True,'major', linewidth=0.5, linestyle='-')

        ax4.grid(True)

        ax5.set_yticks(range(0,max_coord+1,2))
        ax5.set_ylim([-0.5,max_coord+0.5])
        ax5.yaxis.set_major_locator(MultipleLocator(2))
        ax5.yaxis.set_minor_locator(MultipleLocator(1))
        ax5.xaxis.set_major_locator(MultipleLocator(50))
        ax5.xaxis.set_minor_locator(MultipleLocator(25))
        ax5.set_axisbelow(True)
        ax5.yaxis.grid(True,'minor')
        ax5.yaxis.grid(True,'major', linewidth=0.5, linestyle='-')
        ax5.xaxis.grid(True,'minor')
        ax5.xaxis.grid(True,'major', linewidth=0.5, linestyle='-')
        ax5.set_ylabel(coord_cols_map[0]+"\n1st")

        ax6.grid(True)

        ax7.set_yticks(range(0,max_coord+1,2))
        ax7.set_ylim([-0.5,max_coord+0.5])
        ax7.yaxis.set_major_locator(MultipleLocator(2))
        ax7.yaxis.set_minor_locator(MultipleLocator(1))
        ax7.xaxis.set_major_locator(MultipleLocator(50))
        ax7.xaxis.set_minor_locator(MultipleLocator(25))
        ax7.set_axisbelow(True)
        ax7.yaxis.grid(True,'minor')
        ax7.yaxis.grid(True,'major', linewidth=0.5, linestyle='-')
        ax7.xaxis.grid(True,'minor')
        ax7.xaxis.grid(True,'major', linewidth=0.5, linestyle='-')
        ax7.set_ylabel(coord_cols_map[1]+"\n1st")

        ax8.grid(True)

        ax9.set_yticks(range(0,max_coord+1,2))
        ax9.set_ylim([-0.5,max_coord+0.5])
        ax9.yaxis.set_major_locator(MultipleLocator(2))
        ax9.yaxis.set_minor_locator(MultipleLocator(1))
        ax9.xaxis.set_major_locator(MultipleLocator(50))
        ax9.xaxis.set_minor_locator(MultipleLocator(25))
        ax9.set_axisbelow(True)
        ax9.yaxis.grid(True,'minor')
        ax9.yaxis.grid(True,'major', linewidth=0.5, linestyle='-')
        ax9.xaxis.grid(True,'minor')
        ax9.xaxis.grid(True,'major', linewidth=0.5, linestyle='-')
        ax9.set_ylabel(coord_cols_map[0]+"\n2nd")

        ax10.grid(True)

        ax11.set_yticks(range(0,max_coord+1,2))
        ax11.set_ylim([-0.5,max_coord+0.5])
        ax11.yaxis.set_major_locator(MultipleLocator(2))
        ax11.yaxis.set_minor_locator(MultipleLocator(1))
        ax11.xaxis.set_major_locator(MultipleLocator(50))
        ax11.xaxis.set_minor_locator(MultipleLocator(25))
        ax11.set_axisbelow(True)
        ax11.yaxis.grid(True,'minor')
        ax11.yaxis.grid(True,'major', linewidth=0.5, linestyle='-')
        ax11.xaxis.grid(True,'minor')
        ax11.xaxis.grid(True,'major', linewidth=0.5, linestyle='-')
        ax11.set_ylabel(coord_cols_map[1]+"\n2nd")
        ax11.set_xlabel("Time (ns)")

        ax12.grid(True)

        ax13.set_yticks(range(0,max_coord+1,2))
        ax13.set_ylim([-0.5,max_coord+0.5])
        ax13.yaxis.set_major_locator(MultipleLocator(2))
        ax13.yaxis.set_minor_locator(MultipleLocator(1))
        ax13.xaxis.set_major_locator(MultipleLocator(50))
        ax13.xaxis.set_minor_locator(MultipleLocator(25))
        ax13.set_axisbelow(True)
        ax13.yaxis.grid(True,'minor')
        ax13.yaxis.grid(True,'major', linewidth=0.5, linestyle='-')
        ax13.xaxis.grid(True,'minor')
        ax13.xaxis.grid(True,'major', linewidth=0.5, linestyle='-')
        ax13.set_ylabel(coord_cols_map[0]+"\nBoth")

        ax14.grid(True)

        ax15.set_yticks(range(0,max_coord+1,2))
        ax15.set_ylim([-0.5,max_coord+0.5])
        ax15.set_ylabel(coord_cols_map[1]+"\nBoth")
        ax15.yaxis.set_major_locator(MultipleLocator(2))
        ax15.yaxis.set_minor_locator(MultipleLocator(1))
        ax15.xaxis.set_major_locator(MultipleLocator(50))
        ax15.xaxis.set_minor_locator(MultipleLocator(25))
        ax15.set_axisbelow(True)
        ax15.yaxis.grid(True,'minor')
        ax15.yaxis.grid(True,'major', linewidth=0.5, linestyle='-')
        ax15.xaxis.grid(True,'minor')
        ax15.xaxis.grid(True,'major', linewidth=0.5, linestyle='-')

        ax16.grid(True)

        ax17.set_ylabel("R (nm)")
        ax17.yaxis.set_major_locator(MultipleLocator(0.2))
        ax17.yaxis.set_minor_locator(MultipleLocator(0.1))
        ax17.xaxis.set_major_locator(MultipleLocator(50))
        ax17.xaxis.set_minor_locator(MultipleLocator(25))
        ax17.set_axisbelow(True)
        ax17.yaxis.grid(True,'minor')
        ax17.yaxis.grid(True,'major', linewidth=0.5, linestyle='-')
        ax17.xaxis.grid(True,'minor')
        ax17.xaxis.grid(True,'major', linewidth=0.5, linestyle='-')
        ax17.set_ylim([0, 0.6])
        ax17.set_yticks([x/100.0 for x in range(0,50,20)])

        ax17.set_xlabel("Time (ns)")

        ax18.grid(True)
        ax18.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))

        plt.subplots_adjust(hspace = 0.1, wspace = 0.02,
                            left = 0.1, bottom = 0.08,
                            right = 0.95, top = 0.93)
        plt.suptitle(plot_title+" - Traj "+str(traj_id))

        if prefix != None:
            plt.savefig(prefix+"_"+str(traj_id)+".pdf")
        else:
            plt.show()

        plt.close()

    return True

# This plots a stacked timeseries of trajectory properties in the main column
# and histograms of this data in a second column. Sorry for the god-awful
# number of arguments...
def plot_sf_w_3res_timeseries(channel_occ, channel_counts,
                              ion_ts_1st, ion_histo_1st,
                              ion_ts_2nd, ion_histo_2nd,
                              ion_ts_both, ion_histo_both,
                              ion_ts_d1_1st, ion_histo_d1_1st,
                              ion_ts_d2_1st, ion_histo_d2_1st,
                              ion_ts_d3_1st, ion_histo_d3_1st,
                              ion_ts_d1_2nd, ion_histo_d1_2nd,
                              ion_ts_d2_2nd, ion_histo_d2_2nd,
                              ion_ts_d3_2nd, ion_histo_d3_2nd,
                              ion_ts_d1_both, ion_histo_d1_both,
                              ion_ts_d2_both, ion_histo_d2_both,
                              ion_ts_d3_both, ion_histo_d3_both,
                              ion_ts_r, ion_histo_r,
                              time_conv=0.02,
                              prefix=None,
                              plot_title="Stacked Timeseries",
                              max_coord=5,
                              max_length=500,
                              data_skip=10,
                              hist_scale=0.4):

    # This iterates over all the trajectory id's that you computed data for.
    # ion_timeseries[0] is the time values array, but any index would suffice.
    for traj_id in ion_ts_1st[0].keys():

        fig = plt.figure()

        # The grid is a 2 column figure, with 2 + 6 rows.
        gs = gridspec.GridSpec(11, 2,
                               width_ratios=[4,1],
                               height_ratios=[4,1,1,1,1,1,1,1,1,1,2])
        mpl.rcParams['axes.color_cycle'] = ['r', 'g', 'b', 'm', 'c', 'y', 'k']

        # Here we initialize all our subplots on the grid defined above.
        ax1 = plt.subplot(gs[0])
        ax2 = plt.subplot(gs[1], sharey = ax1)
        ax3 = plt.subplot(gs[2], sharex = ax1)
        ax4 = plt.subplot(gs[3], sharey = ax3, sharex = ax2)
        ax5 = plt.subplot(gs[4], sharex = ax1)
        ax6 = plt.subplot(gs[5], sharey = ax3, sharex = ax2)
        ax7 = plt.subplot(gs[6], sharex = ax1)
        ax8 = plt.subplot(gs[7], sharey = ax3, sharex = ax2)
        ax9 = plt.subplot(gs[8], sharex = ax1)
        ax10 = plt.subplot(gs[9], sharey = ax3, sharex = ax2)
        ax11 = plt.subplot(gs[10], sharex = ax1)
        ax12 = plt.subplot(gs[11], sharey = ax3, sharex = ax2)
        ax13 = plt.subplot(gs[12], sharex = ax1)
        ax14 = plt.subplot(gs[13], sharey = ax3, sharex = ax2)
        ax15 = plt.subplot(gs[14], sharex = ax1)
        ax16 = plt.subplot(gs[15], sharey = ax3, sharex = ax2)
        ax17 = plt.subplot(gs[16], sharex = ax1)
        ax18 = plt.subplot(gs[17], sharey = ax3, sharex = ax2)
        ax19 = plt.subplot(gs[18], sharex = ax1)
        ax20 = plt.subplot(gs[19], sharey = ax3, sharex = ax2)
        ax21 = plt.subplot(gs[20], sharex = ax1)
        ax22 = plt.subplot(gs[21], sharey = ax21, sharex = ax2)

        # Many axis are removed.
        plt.setp(ax1.get_xticklabels(), visible=False)
        plt.setp(ax2.get_xticklabels(), visible=False)
        plt.setp(ax2.get_yticklabels(), visible=False)
        plt.setp(ax3.get_xticklabels(), visible=False)
        plt.setp(ax4.get_yticklabels(), visible=False)
        plt.setp(ax4.get_xticklabels(), visible=False)
        plt.setp(ax5.get_xticklabels(), visible=False)
        plt.setp(ax6.get_yticklabels(), visible=False)
        plt.setp(ax6.get_xticklabels(), visible=False)
        plt.setp(ax7.get_xticklabels(), visible=False)
        plt.setp(ax8.get_yticklabels(), visible=False)
        plt.setp(ax8.get_xticklabels(), visible=False)
        plt.setp(ax9.get_xticklabels(), visible=False)
        plt.setp(ax10.get_yticklabels(), visible=False)
        plt.setp(ax10.get_xticklabels(), visible=False)
        plt.setp(ax11.get_xticklabels(), visible=False)
        plt.setp(ax12.get_yticklabels(), visible=False)
        plt.setp(ax12.get_xticklabels(), visible=False)
        plt.setp(ax13.get_xticklabels(), visible=False)
        plt.setp(ax14.get_yticklabels(), visible=False)
        plt.setp(ax14.get_xticklabels(), visible=False)
        plt.setp(ax15.get_xticklabels(), visible=False)
        plt.setp(ax16.get_yticklabels(), visible=False)
        plt.setp(ax16.get_xticklabels(), visible=False)
        plt.setp(ax18.get_xticklabels(), visible=False)
        plt.setp(ax18.get_yticklabels(), visible=False)
        plt.setp(ax19.get_xticklabels(), visible=False)
        plt.setp(ax20.get_xticklabels(), visible=False)
        plt.setp(ax20.get_yticklabels(), visible=False)
        plt.setp(ax22.get_xticklabels(), visible=False)
        plt.setp(ax22.get_yticklabels(), visible=False)

        # Plot 1 - Ion Timeseries
        colors = ['r', 'g', 'b', 'm', 'c', 'y', 'k']
        for color_id, ion_id in enumerate(sorted(ion_ts_1st[0][traj_id].keys(),
                                                 reverse=False)):
            ion_ts_x = ion_ts_1st[1][traj_id][ion_id][::int(data_skip/2)]
            ion_ts_y = ion_ts_1st[0][traj_id][ion_id][::int(data_skip/2)]
            ion_ts_x_scale = [time_conv*pos for pos in ion_ts_x]
            ion_ts_y_flip = [-1*pos for pos in ion_ts_y]
            ax1.scatter(ion_ts_x_scale, ion_ts_y_flip, s=1, color=colors[color_id])

        # Plot 2 - Ion Position Histogram
        for vals, ion_z in zip(ion_histo_1st[0][traj_id], ion_histo_1st[1][traj_id]):
            mpl.rcParams['axes.color_cycle'] = ['r', 'g', 'b', 'm', 'c', 'y', 'k']
            diff = (ion_z[1]-ion_z[0])/2.0
            ion_z_shifted = [-1*(ion-diff) for ion in ion_z[1:]]
            vals_scaled = [val/10.0 for val in vals]
            ax2.plot(vals_scaled, ion_z_shifted)

        # Plot 3 - D1 Coordination Timeseries Per Atom
        for ion_id in sorted(ion_ts_1st[0][traj_id].keys(), reverse=False):
            mpl.rcParams['axes.color_cycle'] = ['r', 'g', 'b', 'm', 'c', 'y', 'k']
            ion_ts_x = ion_ts_d1_1st[1][traj_id][ion_id][::data_skip]
            ion_ts_y = ion_ts_d1_1st[0][traj_id][ion_id][::data_skip]
            ion_ts_x_scale = [time_conv*pos for pos in ion_ts_x]
            ion_ts_y_flip = [1*pos for pos in ion_ts_y]
            ax3.plot(ion_ts_x_scale, ion_ts_y_flip)

        # Plot 4 - Coordination Histogram
        #mpl.rcParams['axes.color_cycle'] = ['r', 'g', 'b', 'm', 'c', 'y', 'k']
        left = 0
        for color_id, hist_data in enumerate(zip(ion_histo_d1_1st[0][traj_id],
                                                      ion_histo_d1_1st[1][traj_id])):
            vals_scaled = array([val*hist_scale for val in hist_data[0]])
            labels_scaled = array(hist_data[1][1:])-1
            ax4.barh(labels_scaled, vals_scaled, align='center',
                     color=colors[color_id], left=left)
            left += vals_scaled

        # Plot 5 - D2 Coordination Timeseries Per Atom
        for ion_id in sorted(ion_ts_1st[0][traj_id].keys(), reverse=False):
            mpl.rcParams['axes.color_cycle'] = ['r', 'g', 'b', 'm', 'c', 'y', 'k']
            ion_ts_x = ion_ts_d1_2nd[1][traj_id][ion_id][::data_skip]
            ion_ts_y = ion_ts_d1_2nd[0][traj_id][ion_id][::data_skip]
            ion_ts_x_scale = [time_conv*pos for pos in ion_ts_x]
            ion_ts_y_flip = [1*pos for pos in ion_ts_y]
            ax5.plot(ion_ts_x_scale, ion_ts_y_flip)

        # Plot 6 - Coordination Histogram
        left = 0
        for color_id, hist_data in enumerate(zip(ion_histo_d1_2nd[0][traj_id],
                                                      ion_histo_d1_2nd[1][traj_id])):
            vals_scaled = array([val*hist_scale for val in hist_data[0]])
            labels_scaled = array(hist_data[1][1:])-1
            ax6.barh(labels_scaled, vals_scaled, align='center',
                     color=colors[color_id], left=left)
            left += vals_scaled

        # Plot 7 - D3 Coordination Timeseries Per Atom
        for ion_id in sorted(ion_ts_1st[0][traj_id].keys(), reverse=False):
            mpl.rcParams['axes.color_cycle'] = ['r', 'g', 'b', 'm', 'c', 'y', 'k']
            ion_ts_x = ion_ts_d1_both[1][traj_id][ion_id][::data_skip]
            ion_ts_y = ion_ts_d1_both[0][traj_id][ion_id][::data_skip]
            ion_ts_x_scale = [time_conv*pos for pos in ion_ts_x]
            ion_ts_y_flip = [1*pos for pos in ion_ts_y]
            ax7.plot(ion_ts_x_scale, ion_ts_y_flip)

        # Plot 8 - Coordination Histogram
        left = 0
        for color_id, hist_data in enumerate(zip(ion_histo_d1_both[0][traj_id],
                                                      ion_histo_d1_both[1][traj_id])):
            vals_scaled = array([val*hist_scale for val in hist_data[0]])
            labels_scaled = array(hist_data[1][1:])-1
            ax8.barh(labels_scaled, vals_scaled, align='center',
                     color=colors[color_id], left=left)
            left += vals_scaled

        # Plot 9 - D1 Coordination Timeseries Per Atom
        for ion_id in sorted(ion_ts_2nd[0][traj_id].keys(), reverse=False):
            mpl.rcParams['axes.color_cycle'] = ['r', 'g', 'b', 'm', 'c', 'y', 'k']
            ion_ts_x = ion_ts_d2_1st[1][traj_id][ion_id][::data_skip]
            ion_ts_y = ion_ts_d2_1st[0][traj_id][ion_id][::data_skip]
            ion_ts_x_scale = [time_conv*pos for pos in ion_ts_x]
            ion_ts_y_flip = [1*pos for pos in ion_ts_y]
            ax9.plot(ion_ts_x_scale, ion_ts_y_flip)

        # Plot 10 - Coordination Histogram
        left = 0
        for color_id, hist_data in enumerate(zip(ion_histo_d2_1st[0][traj_id],
                                                      ion_histo_d2_1st[1][traj_id])):
            vals_scaled = array([val*hist_scale for val in hist_data[0]])
            labels_scaled = array(hist_data[1][1:])-1
            ax10.barh(labels_scaled, vals_scaled, align='center',
                     color=colors[color_id], left=left)
            left += vals_scaled

        # Plot 11 - D2 Coordination Timeseries Per Atom
        for ion_id in sorted(ion_ts_2nd[0][traj_id].keys(), reverse=False):
            mpl.rcParams['axes.color_cycle'] = ['r', 'g', 'b', 'm', 'c', 'y', 'k']
            ion_ts_x = ion_ts_d2_2nd[1][traj_id][ion_id][::data_skip]
            ion_ts_y = ion_ts_d2_2nd[0][traj_id][ion_id][::data_skip]
            ion_ts_x_scale = [time_conv*pos for pos in ion_ts_x]
            ion_ts_y_flip = [1*pos for pos in ion_ts_y]
            ax11.plot(ion_ts_x_scale, ion_ts_y_flip)

        # Plot 12 - Coordination Histogram
        left = 0
        for color_id, hist_data in enumerate(zip(ion_histo_d2_2nd[0][traj_id],
                                                      ion_histo_d2_2nd[1][traj_id])):
            vals_scaled = array([val*hist_scale for val in hist_data[0]])
            labels_scaled = array(hist_data[1][1:])-1
            ax12.barh(labels_scaled, vals_scaled, align='center',
                     color=colors[color_id], left=left)
            left += vals_scaled

        # Plot 13 - D3 Coordination Timeseries Per Atom
        for ion_id in sorted(ion_ts_both[0][traj_id].keys(), reverse=False):
            mpl.rcParams['axes.color_cycle'] = ['r', 'g', 'b', 'm', 'c', 'y', 'k']
            ion_ts_x = ion_ts_d2_both[1][traj_id][ion_id][::data_skip]
            ion_ts_y = ion_ts_d2_both[0][traj_id][ion_id][::data_skip]
            ion_ts_x_scale = [time_conv*pos for pos in ion_ts_x]
            ion_ts_y_flip = [1*pos for pos in ion_ts_y]
            ax13.plot(ion_ts_x_scale, ion_ts_y_flip)

        # Plot 14 - Coordination Histogram
        left = 0
        for color_id, hist_data in enumerate(zip(ion_histo_d2_both[0][traj_id],
                                                      ion_histo_d2_both[1][traj_id])):
            vals_scaled = array([val*hist_scale for val in hist_data[0]])
            labels_scaled = array(hist_data[1][1:])-1
            ax14.barh(labels_scaled, vals_scaled, align='center',
                     color=colors[color_id], left=left)
            left += vals_scaled

        # Plot 15 - D1 Coordination Timeseries Per Atom
        for ion_id in sorted(ion_ts_both[0][traj_id].keys(), reverse=False):
            mpl.rcParams['axes.color_cycle'] = ['r', 'g', 'b', 'm', 'c', 'y', 'k']
            ion_ts_x = ion_ts_d3_1st[1][traj_id][ion_id][::data_skip]
            ion_ts_y = ion_ts_d3_1st[0][traj_id][ion_id][::data_skip]
            ion_ts_x_scale = [time_conv*pos for pos in ion_ts_x]
            ion_ts_y_flip = [1*pos for pos in ion_ts_y]
            ax15.plot(ion_ts_x_scale, ion_ts_y_flip)

        # Plot 16 - Coordination Histogram
        left = 0
        for color_id, hist_data in enumerate(zip(ion_histo_d3_1st[0][traj_id],
                                                      ion_histo_d3_1st[1][traj_id])):
            vals_scaled = array([val*hist_scale for val in hist_data[0]])
            labels_scaled = array(hist_data[1][1:])-1
            ax16.barh(labels_scaled, vals_scaled, align='center',
                     color=colors[color_id], left=left)
            left += vals_scaled

        # Plot 17 - D2 Coordination Timeseries Per Atom
        for ion_id in sorted(ion_ts_both[0][traj_id].keys(), reverse=False):
            mpl.rcParams['axes.color_cycle'] = ['r', 'g', 'b', 'm', 'c', 'y', 'k']
            ion_ts_x = ion_ts_d3_2nd[1][traj_id][ion_id][::data_skip]
            ion_ts_y = ion_ts_d3_2nd[0][traj_id][ion_id][::data_skip]
            ion_ts_x_scale = [time_conv*pos for pos in ion_ts_x]
            ion_ts_y_flip = [1*pos for pos in ion_ts_y]
            ax17.plot(ion_ts_x_scale, ion_ts_y_flip)

        # Plot 18 - Coordination Histogram
        left = 0
        for color_id, hist_data in enumerate(zip(ion_histo_d3_2nd[0][traj_id],
                                                      ion_histo_d3_2nd[1][traj_id])):
            vals_scaled = array([val*hist_scale for val in hist_data[0]])
            labels_scaled = array(hist_data[1][1:])-1
            ax18.barh(labels_scaled, vals_scaled, align='center',
                     color=colors[color_id], left=left)
            left += vals_scaled

        # Plot 19 - D3 Coordination Timeseries Per Atom
        for ion_id in sorted(ion_ts_both[0][traj_id].keys(), reverse=False):
            mpl.rcParams['axes.color_cycle'] = ['r', 'g', 'b', 'm', 'c', 'y', 'k']
            ion_ts_x = ion_ts_d3_both[1][traj_id][ion_id][::data_skip]
            ion_ts_y = ion_ts_d3_both[0][traj_id][ion_id][::data_skip]
            ion_ts_x_scale = [time_conv*pos for pos in ion_ts_x]
            ion_ts_y_flip = [1*pos for pos in ion_ts_y]
            ax19.plot(ion_ts_x_scale, ion_ts_y_flip)

        # Plot 20 - Coordination Histogram
        left = 0
        for color_id, hist_data in enumerate(zip(ion_histo_d3_both[0][traj_id],
                                                      ion_histo_d3_both[1][traj_id])):
            vals_scaled = array([val*hist_scale for val in hist_data[0]])
            labels_scaled = array(hist_data[1][1:])-1
            ax20.barh(labels_scaled, vals_scaled, align='center',
                     color=colors[color_id], left=left)
            left += vals_scaled

        # Plot 21 - Ion R Timeseries
        colors = ['r', 'g', 'b', 'm', 'c', 'y', 'k']
        for color_id, ion_id in enumerate(sorted(ion_ts_r[0][traj_id].keys(),
                                                 reverse=False)):
            ion_ts_x = ion_ts_r[1][traj_id][ion_id][::data_skip/2]
            ion_ts_y = ion_ts_r[0][traj_id][ion_id][::data_skip/2]
            ion_ts_x_scale = [time_conv*pos for pos in ion_ts_x]
            ion_ts_y_flip = [1*pos for pos in ion_ts_y]
            ax21.scatter(ion_ts_x_scale, ion_ts_y_flip, s=1, color=colors[color_id])

        # Plot 22 - Ion R Position Histogram
        for vals, ion_z in zip(ion_histo_r[0][traj_id], ion_histo_r[1][traj_id]):
            mpl.rcParams['axes.color_cycle'] = ['r', 'g', 'b', 'm', 'c', 'y', 'k']
            diff = (ion_z[1]-ion_z[0])/2.0
            ion_z_shifted = [1*(ion-diff) for ion in ion_z[1:]]
            vals_scaled = [val/20.0 for val in vals]
            ax22.plot(vals_scaled, ion_z_shifted)

        ax1.set_ylabel("Axial pos / nm")
        ax1.set_ylim([-0.5, 1.2])
        ax1.set_yticks([x/10.0 for x in range(-4,11,4)])
        ax1.set_xlim([-5,max_length+5])
        ax1.set_xticks(range(0, int(50*round(max_length/50)), 50))
        ax1.yaxis.set_major_locator(MultipleLocator(0.2))
        ax1.yaxis.set_minor_locator(MultipleLocator(0.1))
        ax1.xaxis.set_major_locator(MultipleLocator(50))
        ax1.xaxis.set_minor_locator(MultipleLocator(25))
        ax1.set_axisbelow(True)
        ax1.yaxis.grid(True,'minor')
        ax1.yaxis.grid(True,'major', linewidth=0.5, linestyle='-')
        ax1.xaxis.grid(True,'minor')
        ax1.xaxis.grid(True,'major', linewidth=0.5, linestyle='-')

        ax2.grid(True)
        ax2.set_xlim([0,1.0])
        ax2.set_xticks([x/10.0 for x in range(0,11,2)])

        ax3.set_ylabel("1st")
        ax3.set_yticks(range(0,max_coord+1,2))
        ax3.set_ylim([-0.5,max_coord+0.5])
        ax3.yaxis.set_major_locator(MultipleLocator(2))
        ax3.yaxis.set_minor_locator(MultipleLocator(1))
        ax3.xaxis.set_major_locator(MultipleLocator(50))
        ax3.xaxis.set_minor_locator(MultipleLocator(25))
        ax3.set_axisbelow(True)
        ax3.yaxis.grid(True,'minor')
        ax3.yaxis.grid(True,'major', linewidth=0.5, linestyle='-')
        ax3.xaxis.grid(True,'minor')
        ax3.xaxis.grid(True,'major', linewidth=0.5, linestyle='-')

        ax4.grid(True)

        ax5.set_yticks(range(0,max_coord+1,2))
        ax5.set_ylim([-0.5,max_coord+0.5])
        ax5.yaxis.set_major_locator(MultipleLocator(2))
        ax5.yaxis.set_minor_locator(MultipleLocator(1))
        ax5.xaxis.set_major_locator(MultipleLocator(50))
        ax5.xaxis.set_minor_locator(MultipleLocator(25))
        ax5.set_axisbelow(True)
        ax5.yaxis.grid(True,'minor')
        ax5.yaxis.grid(True,'major', linewidth=0.5, linestyle='-')
        ax5.xaxis.grid(True,'minor')
        ax5.xaxis.grid(True,'major', linewidth=0.5, linestyle='-')
        ax5.set_ylabel("181\n2nd")

        ax6.grid(True)

        ax7.set_yticks(range(0,max_coord+1,2))
        ax7.set_ylim([-0.5,max_coord+0.5])
        ax7.yaxis.set_major_locator(MultipleLocator(2))
        ax7.yaxis.set_minor_locator(MultipleLocator(1))
        ax7.xaxis.set_major_locator(MultipleLocator(50))
        ax7.xaxis.set_minor_locator(MultipleLocator(25))
        ax7.set_axisbelow(True)
        ax7.yaxis.grid(True,'minor')
        ax7.yaxis.grid(True,'major', linewidth=0.5, linestyle='-')
        ax7.xaxis.grid(True,'minor')
        ax7.xaxis.grid(True,'major', linewidth=0.5, linestyle='-')
        ax7.set_ylabel("Bth")

        ax8.grid(True)

        ax9.set_yticks(range(0,max_coord+1,2))
        ax9.set_ylim([-0.5,max_coord+0.5])
        ax9.yaxis.set_major_locator(MultipleLocator(2))
        ax9.yaxis.set_minor_locator(MultipleLocator(1))
        ax9.xaxis.set_major_locator(MultipleLocator(50))
        ax9.xaxis.set_minor_locator(MultipleLocator(25))
        ax9.set_axisbelow(True)
        ax9.yaxis.grid(True,'minor')
        ax9.yaxis.grid(True,'major', linewidth=0.5, linestyle='-')
        ax9.xaxis.grid(True,'minor')
        ax9.xaxis.grid(True,'major', linewidth=0.5, linestyle='-')
        ax9.set_ylabel("1st")

        ax10.grid(True)

        ax11.set_yticks(range(0,max_coord+1,2))
        ax11.set_ylim([-0.5,max_coord+0.5])
        ax11.yaxis.set_major_locator(MultipleLocator(2))
        ax11.yaxis.set_minor_locator(MultipleLocator(1))
        ax11.xaxis.set_major_locator(MultipleLocator(50))
        ax11.xaxis.set_minor_locator(MultipleLocator(25))
        ax11.set_axisbelow(True)
        ax11.yaxis.grid(True,'minor')
        ax11.yaxis.grid(True,'major', linewidth=0.5, linestyle='-')
        ax11.xaxis.grid(True,'minor')
        ax11.xaxis.grid(True,'major', linewidth=0.5, linestyle='-')
        ax11.set_ylabel("178\n2nd")

        ax12.grid(True)

        ax13.set_yticks(range(0,max_coord+1,2))
        ax13.set_ylim([-0.5,max_coord+0.5])
        ax13.yaxis.set_major_locator(MultipleLocator(2))
        ax13.yaxis.set_minor_locator(MultipleLocator(1))
        ax13.xaxis.set_major_locator(MultipleLocator(50))
        ax13.xaxis.set_minor_locator(MultipleLocator(25))
        ax13.set_axisbelow(True)
        ax13.yaxis.grid(True,'minor')
        ax13.yaxis.grid(True,'major', linewidth=0.5, linestyle='-')
        ax13.xaxis.grid(True,'minor')
        ax13.xaxis.grid(True,'major', linewidth=0.5, linestyle='-')
        ax13.set_ylabel("Bth")

        ax14.grid(True)

        ax15.set_yticks(range(0,max_coord+1,2))
        ax15.set_ylim([-0.5,max_coord+0.5])
        ax15.set_ylabel("1st")
        ax15.yaxis.set_major_locator(MultipleLocator(2))
        ax15.yaxis.set_minor_locator(MultipleLocator(1))
        ax15.xaxis.set_major_locator(MultipleLocator(50))
        ax15.xaxis.set_minor_locator(MultipleLocator(25))
        ax15.set_axisbelow(True)
        ax15.yaxis.grid(True,'minor')
        ax15.yaxis.grid(True,'major', linewidth=0.5, linestyle='-')
        ax15.xaxis.grid(True,'minor')
        ax15.xaxis.grid(True,'major', linewidth=0.5, linestyle='-')

        ax16.grid(True)

        ax17.set_yticks(range(0,max_coord+1,2))
        ax17.set_ylim([-0.5,max_coord+0.5])
        ax17.set_ylabel("177\n2nd")
        ax17.yaxis.set_major_locator(MultipleLocator(2))
        ax17.yaxis.set_minor_locator(MultipleLocator(1))
        ax17.xaxis.set_major_locator(MultipleLocator(50))
        ax17.xaxis.set_minor_locator(MultipleLocator(25))
        ax17.set_axisbelow(True)
        ax17.yaxis.grid(True,'minor')
        ax17.yaxis.grid(True,'major', linewidth=0.5, linestyle='-')
        ax17.xaxis.grid(True,'minor')
        ax17.xaxis.grid(True,'major', linewidth=0.5, linestyle='-')

        ax18.grid(True)

        ax19.set_yticks(range(0,max_coord+1,2))
        ax19.set_ylim([-0.5,max_coord+0.5])
        ax19.set_ylabel("Bth")
        ax19.yaxis.set_major_locator(MultipleLocator(2))
        ax19.yaxis.set_minor_locator(MultipleLocator(1))
        ax19.xaxis.set_major_locator(MultipleLocator(50))
        ax19.xaxis.set_minor_locator(MultipleLocator(25))
        ax19.set_axisbelow(True)
        ax19.yaxis.grid(True,'minor')
        ax19.yaxis.grid(True,'major', linewidth=0.5, linestyle='-')
        ax19.xaxis.grid(True,'minor')
        ax19.xaxis.grid(True,'major', linewidth=0.5, linestyle='-')

        ax20.grid(True)

        ax21.set_ylabel("R / nm")
        ax21.yaxis.set_major_locator(MultipleLocator(0.2))
        ax21.yaxis.set_minor_locator(MultipleLocator(0.1))
        ax21.xaxis.set_major_locator(MultipleLocator(50))
        ax21.xaxis.set_minor_locator(MultipleLocator(25))
        ax21.set_axisbelow(True)
        ax21.yaxis.grid(True,'minor')
        ax21.yaxis.grid(True,'major', linewidth=0.5, linestyle='-')
        ax21.xaxis.grid(True,'minor')
        ax21.xaxis.grid(True,'major', linewidth=0.5, linestyle='-')
        ax21.set_ylim([0, 0.6])
        ax21.set_yticks([x/100.0 for x in range(0,50,20)])

        ax21.set_xlabel("Time / ns")

        ax22.grid(True)
        ax22.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))

        plt.subplots_adjust(hspace = 0.1, wspace = 0.02,
                            left = 0.1, bottom = 0.08,
                            right = 0.95, top = 0.93)
        plt.suptitle(plot_title+" - Traj "+str(traj_id))

        if prefix != None:
            plt.savefig(prefix+"_"+str(traj_id)+".pdf")
        else:
            plt.show()

        plt.close()

    return True

# This plots a stacked timeseries of trajectory properties in the main column
# and histograms of this data in a second column. Sorry for the god-awful
# number of arguments...
def plot_sf_w_4res_timeseries(channel_occ, channel_counts,
                              ion_ts_1st, ion_histo_1st,
                              ion_ts_2nd, ion_histo_2nd,
                              ion_ts_both, ion_histo_both,
                              ion_ts_d1_1st, ion_histo_d1_1st,
                              ion_ts_d2_1st, ion_histo_d2_1st,
                              ion_ts_d3_1st, ion_histo_d3_1st,
                              ion_ts_l_1st, ion_histo_l_1st,
                              ion_ts_d1_2nd, ion_histo_d1_2nd,
                              ion_ts_d2_2nd, ion_histo_d2_2nd,
                              ion_ts_d3_2nd, ion_histo_d3_2nd,
                              ion_ts_l_2nd, ion_histo_l_2nd,
                              ion_ts_d1_both, ion_histo_d1_both,
                              ion_ts_d2_both, ion_histo_d2_both,
                              ion_ts_d3_both, ion_histo_d3_both,
                              ion_ts_l_both, ion_histo_l_both,
                              time_conv=0.02,
                              prefix=None,
                              plot_title="Stacked Timeseries",
                              max_coord=5,
                              max_length=500,
                              data_skip=10,
                              hist_scale=0.4):

    # This iterates over all the trajectory id's that you computed data for.
    # ion_timeseries[0] is the time values array, but any index would suffice.
    for traj_id in ion_ts_1st[0].keys():

        fig = plt.figure()

        # The grid is a 2 column figure, with 2 + 6 rows.
        gs = gridspec.GridSpec(13, 2,
                               width_ratios=[4,1],
                               height_ratios=[4,1,1,1,1,1,1,1,1,1,1,1,1])
        mpl.rcParams['axes.color_cycle'] = ['r', 'g', 'b', 'm', 'c', 'y', 'k']

        # Here we initialize all our subplots on the grid defined above.
        ax1 = plt.subplot(gs[0])
        ax2 = plt.subplot(gs[1], sharey = ax1)
        ax3 = plt.subplot(gs[2], sharex = ax1)
        ax4 = plt.subplot(gs[3], sharey = ax3, sharex = ax2)
        ax5 = plt.subplot(gs[4], sharex = ax1)
        ax6 = plt.subplot(gs[5], sharey = ax3, sharex = ax2)
        ax7 = plt.subplot(gs[6], sharex = ax1)
        ax8 = plt.subplot(gs[7], sharey = ax3, sharex = ax2)
        ax9 = plt.subplot(gs[8], sharex = ax1)
        ax10 = plt.subplot(gs[9], sharey = ax3, sharex = ax2)
        ax11 = plt.subplot(gs[10], sharex = ax1)
        ax12 = plt.subplot(gs[11], sharey = ax3, sharex = ax2)
        ax13 = plt.subplot(gs[12], sharex = ax1)
        ax14 = plt.subplot(gs[13], sharey = ax3, sharex = ax2)
        ax15 = plt.subplot(gs[14], sharex = ax1)
        ax16 = plt.subplot(gs[15], sharey = ax3, sharex = ax2)
        ax17 = plt.subplot(gs[16], sharex = ax1)
        ax18 = plt.subplot(gs[17], sharey = ax3, sharex = ax2)
        ax19 = plt.subplot(gs[18], sharex = ax1)
        ax20 = plt.subplot(gs[19], sharey = ax3, sharex = ax2)
        ax21 = plt.subplot(gs[20], sharex = ax1)
        ax22 = plt.subplot(gs[21], sharey = ax3, sharex = ax2)
        ax23 = plt.subplot(gs[22], sharex = ax1)
        ax24 = plt.subplot(gs[23], sharey = ax3, sharex = ax2)
        ax25 = plt.subplot(gs[24], sharex = ax1)
        ax26 = plt.subplot(gs[25], sharey = ax3, sharex = ax2)


        # Many axis are removed.
        plt.setp(ax1.get_xticklabels(), visible=False)
        plt.setp(ax2.get_xticklabels(), visible=False)
        plt.setp(ax2.get_yticklabels(), visible=False)
        plt.setp(ax3.get_xticklabels(), visible=False)
        plt.setp(ax4.get_yticklabels(), visible=False)
        plt.setp(ax4.get_xticklabels(), visible=False)
        plt.setp(ax5.get_xticklabels(), visible=False)
        plt.setp(ax6.get_yticklabels(), visible=False)
        plt.setp(ax6.get_xticklabels(), visible=False)
        plt.setp(ax7.get_xticklabels(), visible=False)
        plt.setp(ax8.get_yticklabels(), visible=False)
        plt.setp(ax8.get_xticklabels(), visible=False)
        plt.setp(ax9.get_xticklabels(), visible=False)
        plt.setp(ax10.get_yticklabels(), visible=False)
        plt.setp(ax10.get_xticklabels(), visible=False)
        plt.setp(ax11.get_xticklabels(), visible=False)
        plt.setp(ax12.get_yticklabels(), visible=False)
        plt.setp(ax12.get_xticklabels(), visible=False)
        plt.setp(ax13.get_xticklabels(), visible=False)
        plt.setp(ax14.get_yticklabels(), visible=False)
        plt.setp(ax14.get_xticklabels(), visible=False)
        plt.setp(ax15.get_xticklabels(), visible=False)
        plt.setp(ax16.get_yticklabels(), visible=False)
        plt.setp(ax16.get_xticklabels(), visible=False)
        plt.setp(ax18.get_xticklabels(), visible=False)
        plt.setp(ax18.get_yticklabels(), visible=False)
        plt.setp(ax19.get_xticklabels(), visible=False)
        plt.setp(ax20.get_xticklabels(), visible=False)
        plt.setp(ax20.get_yticklabels(), visible=False)
        plt.setp(ax21.get_xticklabels(), visible=False)
        plt.setp(ax22.get_xticklabels(), visible=False)
        plt.setp(ax22.get_yticklabels(), visible=False)
        plt.setp(ax23.get_xticklabels(), visible=False)
        plt.setp(ax24.get_xticklabels(), visible=False)
        plt.setp(ax24.get_yticklabels(), visible=False)
        plt.setp(ax26.get_xticklabels(), visible=False)
        plt.setp(ax26.get_yticklabels(), visible=False)

        # Plot 1 - Ion Timeseries
        colors = ['r', 'g', 'b', 'm', 'c', 'y', 'k']
        for color_id, ion_id in enumerate(sorted(ion_ts_1st[0][traj_id].keys(),
                                                 reverse=False)):
            ion_ts_x = ion_ts_1st[1][traj_id][ion_id][::int(data_skip/2)]
            ion_ts_y = ion_ts_1st[0][traj_id][ion_id][::int(data_skip/2)]
            ion_ts_x_scale = [time_conv*pos for pos in ion_ts_x]
            ion_ts_y_flip = [-1*pos for pos in ion_ts_y]
            ax1.scatter(ion_ts_x_scale, ion_ts_y_flip, s=1, color=colors[color_id])

        # Plot 2 - Ion Position Histogram
        for vals, ion_z in zip(ion_histo_1st[0][traj_id], ion_histo_1st[1][traj_id]):
            mpl.rcParams['axes.color_cycle'] = ['r', 'g', 'b', 'm', 'c', 'y', 'k']
            diff = (ion_z[1]-ion_z[0])/2.0
            ion_z_shifted = [-1*(ion-diff) for ion in ion_z[1:]]
            vals_scaled = [val/10.0 for val in vals]
            ax2.plot(vals_scaled, ion_z_shifted)

        # Plot 3 - D1 Coordination Timeseries Per Atom
        for ion_id in sorted(ion_ts_1st[0][traj_id].keys(), reverse=False):
            mpl.rcParams['axes.color_cycle'] = ['r', 'g', 'b', 'm', 'c', 'y', 'k']
            ion_ts_x = ion_ts_d1_1st[1][traj_id][ion_id][::data_skip]
            ion_ts_y = ion_ts_d1_1st[0][traj_id][ion_id][::data_skip]
            ion_ts_x_scale = [time_conv*pos for pos in ion_ts_x]
            ion_ts_y_flip = [1*pos for pos in ion_ts_y]
            ax3.plot(ion_ts_x_scale, ion_ts_y_flip)

        # Plot 4 - Coordination Histogram
        #mpl.rcParams['axes.color_cycle'] = ['r', 'g', 'b', 'm', 'c', 'y', 'k']
        left = 0
        for color_id, hist_data in enumerate(zip(ion_histo_d1_1st[0][traj_id],
                                                      ion_histo_d1_1st[1][traj_id])):
            vals_scaled = array([val*hist_scale for val in hist_data[0]])
            labels_scaled = array(hist_data[1][1:])-1
            ax4.barh(labels_scaled, vals_scaled, align='center',
                     color=colors[color_id], left=left)
            left += vals_scaled

        # Plot 5 - D2 Coordination Timeseries Per Atom
        for ion_id in sorted(ion_ts_1st[0][traj_id].keys(), reverse=False):
            mpl.rcParams['axes.color_cycle'] = ['r', 'g', 'b', 'm', 'c', 'y', 'k']
            ion_ts_x = ion_ts_d1_2nd[1][traj_id][ion_id][::data_skip]
            ion_ts_y = ion_ts_d1_2nd[0][traj_id][ion_id][::data_skip]
            ion_ts_x_scale = [time_conv*pos for pos in ion_ts_x]
            ion_ts_y_flip = [1*pos for pos in ion_ts_y]
            ax5.plot(ion_ts_x_scale, ion_ts_y_flip)

        # Plot 6 - Coordination Histogram
        left = 0
        for color_id, hist_data in enumerate(zip(ion_histo_d1_2nd[0][traj_id],
                                                      ion_histo_d1_2nd[1][traj_id])):
            vals_scaled = array([val*hist_scale for val in hist_data[0]])
            labels_scaled = array(hist_data[1][1:])-1
            ax6.barh(labels_scaled, vals_scaled, align='center',
                     color=colors[color_id], left=left)
            left += vals_scaled

        # Plot 7 - D3 Coordination Timeseries Per Atom
        for ion_id in sorted(ion_ts_1st[0][traj_id].keys(), reverse=False):
            mpl.rcParams['axes.color_cycle'] = ['r', 'g', 'b', 'm', 'c', 'y', 'k']
            ion_ts_x = ion_ts_d1_both[1][traj_id][ion_id][::data_skip]
            ion_ts_y = ion_ts_d1_both[0][traj_id][ion_id][::data_skip]
            ion_ts_x_scale = [time_conv*pos for pos in ion_ts_x]
            ion_ts_y_flip = [1*pos for pos in ion_ts_y]
            ax7.plot(ion_ts_x_scale, ion_ts_y_flip)

        # Plot 8 - Coordination Histogram
        left = 0
        for color_id, hist_data in enumerate(zip(ion_histo_d1_both[0][traj_id],
                                                      ion_histo_d1_both[1][traj_id])):
            vals_scaled = array([val*hist_scale for val in hist_data[0]])
            labels_scaled = array(hist_data[1][1:])-1
            ax8.barh(labels_scaled, vals_scaled, align='center',
                     color=colors[color_id], left=left)
            left += vals_scaled

        # Plot 9 - D1 Coordination Timeseries Per Atom
        for ion_id in sorted(ion_ts_2nd[0][traj_id].keys(), reverse=False):
            mpl.rcParams['axes.color_cycle'] = ['r', 'g', 'b', 'm', 'c', 'y', 'k']
            ion_ts_x = ion_ts_d2_1st[1][traj_id][ion_id][::data_skip]
            ion_ts_y = ion_ts_d2_1st[0][traj_id][ion_id][::data_skip]
            ion_ts_x_scale = [time_conv*pos for pos in ion_ts_x]
            ion_ts_y_flip = [1*pos for pos in ion_ts_y]
            ax9.plot(ion_ts_x_scale, ion_ts_y_flip)

        # Plot 10 - Coordination Histogram
        left = 0
        for color_id, hist_data in enumerate(zip(ion_histo_d2_1st[0][traj_id],
                                                      ion_histo_d2_1st[1][traj_id])):
            vals_scaled = array([val*hist_scale for val in hist_data[0]])
            labels_scaled = array(hist_data[1][1:])-1
            ax10.barh(labels_scaled, vals_scaled, align='center',
                     color=colors[color_id], left=left)
            left += vals_scaled

        # Plot 11 - D2 Coordination Timeseries Per Atom
        for ion_id in sorted(ion_ts_2nd[0][traj_id].keys(), reverse=False):
            mpl.rcParams['axes.color_cycle'] = ['r', 'g', 'b', 'm', 'c', 'y', 'k']
            ion_ts_x = ion_ts_d2_2nd[1][traj_id][ion_id][::data_skip]
            ion_ts_y = ion_ts_d2_2nd[0][traj_id][ion_id][::data_skip]
            ion_ts_x_scale = [time_conv*pos for pos in ion_ts_x]
            ion_ts_y_flip = [1*pos for pos in ion_ts_y]
            ax11.plot(ion_ts_x_scale, ion_ts_y_flip)

        # Plot 12 - Coordination Histogram
        left = 0
        for color_id, hist_data in enumerate(zip(ion_histo_d2_2nd[0][traj_id],
                                                      ion_histo_d2_2nd[1][traj_id])):
            vals_scaled = array([val*hist_scale for val in hist_data[0]])
            labels_scaled = array(hist_data[1][1:])-1
            ax12.barh(labels_scaled, vals_scaled, align='center',
                     color=colors[color_id], left=left)
            left += vals_scaled

        # Plot 13 - D3 Coordination Timeseries Per Atom
        for ion_id in sorted(ion_ts_both[0][traj_id].keys(), reverse=False):
            mpl.rcParams['axes.color_cycle'] = ['r', 'g', 'b', 'm', 'c', 'y', 'k']
            ion_ts_x = ion_ts_d2_both[1][traj_id][ion_id][::data_skip]
            ion_ts_y = ion_ts_d2_both[0][traj_id][ion_id][::data_skip]
            ion_ts_x_scale = [time_conv*pos for pos in ion_ts_x]
            ion_ts_y_flip = [1*pos for pos in ion_ts_y]
            ax13.plot(ion_ts_x_scale, ion_ts_y_flip)

        # Plot 14 - Coordination Histogram
        left = 0
        for color_id, hist_data in enumerate(zip(ion_histo_d2_both[0][traj_id],
                                                      ion_histo_d2_both[1][traj_id])):
            vals_scaled = array([val*hist_scale for val in hist_data[0]])
            labels_scaled = array(hist_data[1][1:])-1
            ax14.barh(labels_scaled, vals_scaled, align='center',
                     color=colors[color_id], left=left)
            left += vals_scaled

        # Plot 15 - D1 Coordination Timeseries Per Atom
        for ion_id in sorted(ion_ts_both[0][traj_id].keys(), reverse=False):
            mpl.rcParams['axes.color_cycle'] = ['r', 'g', 'b', 'm', 'c', 'y', 'k']
            ion_ts_x = ion_ts_d3_1st[1][traj_id][ion_id][::data_skip]
            ion_ts_y = ion_ts_d3_1st[0][traj_id][ion_id][::data_skip]
            ion_ts_x_scale = [time_conv*pos for pos in ion_ts_x]
            ion_ts_y_flip = [1*pos for pos in ion_ts_y]
            ax15.plot(ion_ts_x_scale, ion_ts_y_flip)

        # Plot 16 - Coordination Histogram
        left = 0
        for color_id, hist_data in enumerate(zip(ion_histo_d3_1st[0][traj_id],
                                                      ion_histo_d3_1st[1][traj_id])):
            vals_scaled = array([val*hist_scale for val in hist_data[0]])
            labels_scaled = array(hist_data[1][1:])-1
            ax16.barh(labels_scaled, vals_scaled, align='center',
                     color=colors[color_id], left=left)
            left += vals_scaled

        # Plot 17 - D2 Coordination Timeseries Per Atom
        for ion_id in sorted(ion_ts_both[0][traj_id].keys(), reverse=False):
            mpl.rcParams['axes.color_cycle'] = ['r', 'g', 'b', 'm', 'c', 'y', 'k']
            ion_ts_x = ion_ts_d3_2nd[1][traj_id][ion_id][::data_skip]
            ion_ts_y = ion_ts_d3_2nd[0][traj_id][ion_id][::data_skip]
            ion_ts_x_scale = [time_conv*pos for pos in ion_ts_x]
            ion_ts_y_flip = [1*pos for pos in ion_ts_y]
            ax17.plot(ion_ts_x_scale, ion_ts_y_flip)

        # Plot 18 - Coordination Histogram
        left = 0
        for color_id, hist_data in enumerate(zip(ion_histo_d3_2nd[0][traj_id],
                                                      ion_histo_d3_2nd[1][traj_id])):
            vals_scaled = array([val*hist_scale for val in hist_data[0]])
            labels_scaled = array(hist_data[1][1:])-1
            ax18.barh(labels_scaled, vals_scaled, align='center',
                     color=colors[color_id], left=left)
            left += vals_scaled

        # Plot 19 - D3 Coordination Timeseries Per Atom
        for ion_id in sorted(ion_ts_both[0][traj_id].keys(), reverse=False):
            mpl.rcParams['axes.color_cycle'] = ['r', 'g', 'b', 'm', 'c', 'y', 'k']
            ion_ts_x = ion_ts_d3_both[1][traj_id][ion_id][::data_skip]
            ion_ts_y = ion_ts_d3_both[0][traj_id][ion_id][::data_skip]
            ion_ts_x_scale = [time_conv*pos for pos in ion_ts_x]
            ion_ts_y_flip = [1*pos for pos in ion_ts_y]
            ax19.plot(ion_ts_x_scale, ion_ts_y_flip)

        # Plot 20 - Coordination Histogram
        left = 0
        for color_id, hist_data in enumerate(zip(ion_histo_d3_both[0][traj_id],
                                                      ion_histo_d3_both[1][traj_id])):
            vals_scaled = array([val*hist_scale for val in hist_data[0]])
            labels_scaled = array(hist_data[1][1:])-1
            ax20.barh(labels_scaled, vals_scaled, align='center',
                     color=colors[color_id], left=left)
            left += vals_scaled

        # Plot 15 - L Coordination Timeseries Per Atom
        for ion_id in sorted(ion_ts_both[0][traj_id].keys(), reverse=False):
            mpl.rcParams['axes.color_cycle'] = ['r', 'g', 'b', 'm', 'c', 'y', 'k']
            ion_ts_x = ion_ts_l_1st[1][traj_id][ion_id][::data_skip]
            ion_ts_y = ion_ts_l_1st[0][traj_id][ion_id][::data_skip]
            ion_ts_x_scale = [time_conv*pos for pos in ion_ts_x]
            ion_ts_y_flip = [1*pos for pos in ion_ts_y]
            ax21.plot(ion_ts_x_scale, ion_ts_y_flip)

        # Plot 16 - Coordination Histogram
        left = 0
        for color_id, hist_data in enumerate(zip(ion_histo_l_1st[0][traj_id],
                                                      ion_histo_l_1st[1][traj_id])):
            vals_scaled = array([val*hist_scale for val in hist_data[0]])
            labels_scaled = array(hist_data[1][1:])-1
            ax22.barh(labels_scaled, vals_scaled, align='center',
                     color=colors[color_id], left=left)
            left += vals_scaled

        # Plot 17 - L Coordination Timeseries Per Atom
        for ion_id in sorted(ion_ts_both[0][traj_id].keys(), reverse=False):
            mpl.rcParams['axes.color_cycle'] = ['r', 'g', 'b', 'm', 'c', 'y', 'k']
            ion_ts_x = ion_ts_l_2nd[1][traj_id][ion_id][::data_skip]
            ion_ts_y = ion_ts_l_2nd[0][traj_id][ion_id][::data_skip]
            ion_ts_x_scale = [time_conv*pos for pos in ion_ts_x]
            ion_ts_y_flip = [1*pos for pos in ion_ts_y]
            ax23.plot(ion_ts_x_scale, ion_ts_y_flip)

        # Plot 24 - Coordination Histogram
        left = 0
        for color_id, hist_data in enumerate(zip(ion_histo_l_2nd[0][traj_id],
                                                      ion_histo_l_2nd[1][traj_id])):
            vals_scaled = array([val*hist_scale for val in hist_data[0]])
            labels_scaled = array(hist_data[1][1:])-1
            ax24.barh(labels_scaled, vals_scaled, align='center',
                     color=colors[color_id], left=left)
            left += vals_scaled

        # Plot 25 - L Coordination Timeseries Per Atom
        for ion_id in sorted(ion_ts_both[0][traj_id].keys(), reverse=False):
            mpl.rcParams['axes.color_cycle'] = ['r', 'g', 'b', 'm', 'c', 'y', 'k']
            ion_ts_x = ion_ts_l_both[1][traj_id][ion_id][::data_skip]
            ion_ts_y = ion_ts_l_both[0][traj_id][ion_id][::data_skip]
            ion_ts_x_scale = [time_conv*pos for pos in ion_ts_x]
            ion_ts_y_flip = [1*pos for pos in ion_ts_y]
            ax25.plot(ion_ts_x_scale, ion_ts_y_flip)

        # Plot 26 - Coordination Histogram
        left = 0
        for color_id, hist_data in enumerate(zip(ion_histo_l_both[0][traj_id],
                                                      ion_histo_l_both[1][traj_id])):
            vals_scaled = array([val*hist_scale for val in hist_data[0]])
            labels_scaled = array(hist_data[1][1:])-1
            ax26.barh(labels_scaled, vals_scaled, align='center',
                     color=colors[color_id], left=left)
            left += vals_scaled


        ax1.set_ylabel("Axial pos / nm")
        #ax1.set_ylim([-0.5, 1.2])
        #ax1.set_yticks([x/10.0 for x in range(-4,11,4)])
        ax1.set_ylim([-0.8, 1.2])
        ax1.set_yticks([x/10.0 for x in range(-9,11,4)])
        ax1.set_xlim([-5,max_length+5])
        ax1.set_xticks(range(0, int(50*round(max_length/50)), 50))
        ax1.yaxis.set_major_locator(MultipleLocator(0.2))
        ax1.yaxis.set_minor_locator(MultipleLocator(0.1))
        ax1.xaxis.set_major_locator(MultipleLocator(50))
        ax1.xaxis.set_minor_locator(MultipleLocator(25))
        ax1.set_axisbelow(True)
        ax1.yaxis.grid(True,'minor')
        ax1.yaxis.grid(True,'major', linewidth=0.5, linestyle='-')
        ax1.xaxis.grid(True,'minor')
        ax1.xaxis.grid(True,'major', linewidth=0.5, linestyle='-')

        ax2.grid(True)
        ax2.set_xlim([0,1.0])
        ax2.set_xticks([x/10.0 for x in range(0,11,2)])

        ax3.set_ylabel("1st")
        ax3.set_yticks(range(0,max_coord+1,2))
        ax3.set_ylim([-0.5,max_coord+0.5])
        ax3.yaxis.set_major_locator(MultipleLocator(2))
        ax3.yaxis.set_minor_locator(MultipleLocator(1))
        ax3.xaxis.set_major_locator(MultipleLocator(50))
        ax3.xaxis.set_minor_locator(MultipleLocator(25))
        ax3.set_axisbelow(True)
        ax3.yaxis.grid(True,'minor')
        ax3.yaxis.grid(True,'major', linewidth=0.5, linestyle='-')
        ax3.xaxis.grid(True,'minor')
        ax3.xaxis.grid(True,'major', linewidth=0.5, linestyle='-')

        ax4.grid(True)

        ax5.set_yticks(range(0,max_coord+1,2))
        ax5.set_ylim([-0.5,max_coord+0.5])
        ax5.yaxis.set_major_locator(MultipleLocator(2))
        ax5.yaxis.set_minor_locator(MultipleLocator(1))
        ax5.xaxis.set_major_locator(MultipleLocator(50))
        ax5.xaxis.set_minor_locator(MultipleLocator(25))
        ax5.set_axisbelow(True)
        ax5.yaxis.grid(True,'minor')
        ax5.yaxis.grid(True,'major', linewidth=0.5, linestyle='-')
        ax5.xaxis.grid(True,'minor')
        ax5.xaxis.grid(True,'major', linewidth=0.5, linestyle='-')
        ax5.set_ylabel("181\n2nd")

        ax6.grid(True)

        ax7.set_yticks(range(0,max_coord+1,2))
        ax7.set_ylim([-0.5,max_coord+0.5])
        ax7.yaxis.set_major_locator(MultipleLocator(2))
        ax7.yaxis.set_minor_locator(MultipleLocator(1))
        ax7.xaxis.set_major_locator(MultipleLocator(50))
        ax7.xaxis.set_minor_locator(MultipleLocator(25))
        ax7.set_axisbelow(True)
        ax7.yaxis.grid(True,'minor')
        ax7.yaxis.grid(True,'major', linewidth=0.5, linestyle='-')
        ax7.xaxis.grid(True,'minor')
        ax7.xaxis.grid(True,'major', linewidth=0.5, linestyle='-')
        ax7.set_ylabel("Bth")

        ax8.grid(True)

        ax9.set_yticks(range(0,max_coord+1,2))
        ax9.set_ylim([-0.5,max_coord+0.5])
        ax9.yaxis.set_major_locator(MultipleLocator(2))
        ax9.yaxis.set_minor_locator(MultipleLocator(1))
        ax9.xaxis.set_major_locator(MultipleLocator(50))
        ax9.xaxis.set_minor_locator(MultipleLocator(25))
        ax9.set_axisbelow(True)
        ax9.yaxis.grid(True,'minor')
        ax9.yaxis.grid(True,'major', linewidth=0.5, linestyle='-')
        ax9.xaxis.grid(True,'minor')
        ax9.xaxis.grid(True,'major', linewidth=0.5, linestyle='-')
        ax9.set_ylabel("1st")

        ax10.grid(True)

        ax11.set_yticks(range(0,max_coord+1,2))
        ax11.set_ylim([-0.5,max_coord+0.5])
        ax11.yaxis.set_major_locator(MultipleLocator(2))
        ax11.yaxis.set_minor_locator(MultipleLocator(1))
        ax11.xaxis.set_major_locator(MultipleLocator(50))
        ax11.xaxis.set_minor_locator(MultipleLocator(25))
        ax11.set_axisbelow(True)
        ax11.yaxis.grid(True,'minor')
        ax11.yaxis.grid(True,'major', linewidth=0.5, linestyle='-')
        ax11.xaxis.grid(True,'minor')
        ax11.xaxis.grid(True,'major', linewidth=0.5, linestyle='-')
        ax11.set_ylabel("178\n2nd")

        ax12.grid(True)

        ax13.set_yticks(range(0,max_coord+1,2))
        ax13.set_ylim([-0.5,max_coord+0.5])
        ax13.yaxis.set_major_locator(MultipleLocator(2))
        ax13.yaxis.set_minor_locator(MultipleLocator(1))
        ax13.xaxis.set_major_locator(MultipleLocator(50))
        ax13.xaxis.set_minor_locator(MultipleLocator(25))
        ax13.set_axisbelow(True)
        ax13.yaxis.grid(True,'minor')
        ax13.yaxis.grid(True,'major', linewidth=0.5, linestyle='-')
        ax13.xaxis.grid(True,'minor')
        ax13.xaxis.grid(True,'major', linewidth=0.5, linestyle='-')
        ax13.set_ylabel("Bth")

        ax14.grid(True)

        ax15.set_yticks(range(0,max_coord+1,2))
        ax15.set_ylim([-0.5,max_coord+0.5])
        ax15.set_ylabel("1st")
        ax15.yaxis.set_major_locator(MultipleLocator(2))
        ax15.yaxis.set_minor_locator(MultipleLocator(1))
        ax15.xaxis.set_major_locator(MultipleLocator(50))
        ax15.xaxis.set_minor_locator(MultipleLocator(25))
        ax15.set_axisbelow(True)
        ax15.yaxis.grid(True,'minor')
        ax15.yaxis.grid(True,'major', linewidth=0.5, linestyle='-')
        ax15.xaxis.grid(True,'minor')
        ax15.xaxis.grid(True,'major', linewidth=0.5, linestyle='-')

        ax16.grid(True)

        ax17.set_yticks(range(0,max_coord+1,2))
        ax17.set_ylim([-0.5,max_coord+0.5])
        ax17.set_ylabel("177\n2nd")
        ax17.yaxis.set_major_locator(MultipleLocator(2))
        ax17.yaxis.set_minor_locator(MultipleLocator(1))
        ax17.xaxis.set_major_locator(MultipleLocator(50))
        ax17.xaxis.set_minor_locator(MultipleLocator(25))
        ax17.set_axisbelow(True)
        ax17.yaxis.grid(True,'minor')
        ax17.yaxis.grid(True,'major', linewidth=0.5, linestyle='-')
        ax17.xaxis.grid(True,'minor')
        ax17.xaxis.grid(True,'major', linewidth=0.5, linestyle='-')

        ax18.grid(True)

        ax19.set_yticks(range(0,max_coord+1,2))
        ax19.set_ylim([-0.5,max_coord+0.5])
        ax19.set_ylabel("Bth")
        ax19.yaxis.set_major_locator(MultipleLocator(2))
        ax19.yaxis.set_minor_locator(MultipleLocator(1))
        ax19.xaxis.set_major_locator(MultipleLocator(50))
        ax19.xaxis.set_minor_locator(MultipleLocator(25))
        ax19.set_axisbelow(True)
        ax19.yaxis.grid(True,'minor')
        ax19.yaxis.grid(True,'major', linewidth=0.5, linestyle='-')
        ax19.xaxis.grid(True,'minor')
        ax19.xaxis.grid(True,'major', linewidth=0.5, linestyle='-')

        ax20.grid(True)

        ax21.set_yticks(range(0,max_coord+1,2))
        ax21.set_ylim([-0.5,max_coord+0.5])
        ax21.set_ylabel("1st")
        ax21.yaxis.set_major_locator(MultipleLocator(2))
        ax21.yaxis.set_minor_locator(MultipleLocator(1))
        ax21.xaxis.set_major_locator(MultipleLocator(50))
        ax21.xaxis.set_minor_locator(MultipleLocator(25))
        ax21.set_axisbelow(True)
        ax21.yaxis.grid(True,'minor')
        ax21.yaxis.grid(True,'major', linewidth=0.5, linestyle='-')
        ax21.xaxis.grid(True,'minor')
        ax21.xaxis.grid(True,'major', linewidth=0.5, linestyle='-')

        ax22.grid(True)

        ax23.set_yticks(range(0,max_coord+1,2))
        ax23.set_ylim([-0.5,max_coord+0.5])
        ax23.set_ylabel("L\n2nd")
        ax23.yaxis.set_major_locator(MultipleLocator(2))
        ax23.yaxis.set_minor_locator(MultipleLocator(1))
        ax23.xaxis.set_major_locator(MultipleLocator(50))
        ax23.xaxis.set_minor_locator(MultipleLocator(25))
        ax23.set_axisbelow(True)
        ax23.yaxis.grid(True,'minor')
        ax23.yaxis.grid(True,'major', linewidth=0.5, linestyle='-')
        ax23.xaxis.grid(True,'minor')
        ax23.xaxis.grid(True,'major', linewidth=0.5, linestyle='-')

        ax24.grid(True)

        ax25.set_yticks(range(0,max_coord+1,2))
        ax25.set_ylim([-0.5,max_coord+0.5])
        ax25.set_ylabel("Bth")
        ax25.yaxis.set_major_locator(MultipleLocator(2))
        ax25.yaxis.set_minor_locator(MultipleLocator(1))
        ax25.xaxis.set_major_locator(MultipleLocator(50))
        ax25.xaxis.set_minor_locator(MultipleLocator(25))
        ax25.set_axisbelow(True)
        ax25.yaxis.grid(True,'minor')
        ax25.yaxis.grid(True,'major', linewidth=0.5, linestyle='-')
        ax25.xaxis.grid(True,'minor')
        ax25.xaxis.grid(True,'major', linewidth=0.5, linestyle='-')

        ax25.set_xlabel("Time / ns")

        ax26.grid(True)
        ax26.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))

        plt.subplots_adjust(hspace = 0.1, wspace = 0.02,
                            left = 0.1, bottom = 0.08,
                            right = 0.95, top = 0.93)
        plt.suptitle(plot_title+" - Traj "+str(traj_id))

        if prefix != None:
            plt.savefig(prefix+"_"+str(traj_id)+".pdf")
        else:
            plt.show()

        plt.close()

    return True