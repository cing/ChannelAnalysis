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

# This plots a stacked timeseries of trajectory properties in the main column
# and histograms of this data in a second column. Sorry for the god-awful
# number of arguments...
def plot_navab_figure1(channel_occ, channel_counts,
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
                       plot_title="Figure 1 - Stacked Timeseries",
                       max_coord=4,
                       max_length=500,
                       data_skip=10):

    # This iterates over all the trajectory id's that you computed data for.
    # ion_timeseries[0] is the time values array, but any index would suffice.
    for traj_id in ion_ts_1st[0].keys():

        # The grid is a 2 column figure, with 2 + 6 rows.
        gs = gridspec.GridSpec(9, 2,
                               width_ratios=[4,1],
                               height_ratios=[4,1,1,1,1,1,1,1,1])
        mpl.rcParams['axes.color_cycle'] = ['r', 'g', 'b', 'm', 'c']

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
        colors = ['r', 'g', 'b', 'm', 'c']
        for color_id, ion_id in enumerate(sorted(ion_ts_1st[0][traj_id].keys(),
                                                 reverse=False)):
            ion_ts_x = ion_ts_1st[1][traj_id][ion_id][::int(data_skip/2)]
            ion_ts_y = ion_ts_1st[0][traj_id][ion_id][::int(data_skip/2)]
            ion_ts_x_scale = [time_conv*pos for pos in ion_ts_x]
            ion_ts_y_flip = [-1*pos for pos in ion_ts_y]
            ax1.scatter(ion_ts_x_scale, ion_ts_y_flip, s=1, color=colors[color_id])

        # Plot 2 - Ion Position Histogram
        for vals, ion_z in zip(ion_histo_1st[0][traj_id], ion_histo_1st[1][traj_id]):
            mpl.rcParams['axes.color_cycle'] = ['r', 'g', 'b', 'm', 'c']
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
            mpl.rcParams['axes.color_cycle'] = ['r', 'g', 'b', 'm', 'c']
            ion_ts_x = ion_ts_e_1st[1][traj_id][ion_id][::data_skip]
            ion_ts_y = ion_ts_e_1st[0][traj_id][ion_id][::data_skip]
            ion_ts_x_scale = [time_conv*pos for pos in ion_ts_x]
            ion_ts_y_flip = [1*pos for pos in ion_ts_y]
            ax5.plot(ion_ts_x_scale, ion_ts_y_flip)

        # Plot 6 - Coordination Histogram
        channel_occ_hist_x = ion_histo_e_1st[1][traj_id]
        channel_occ_hist_y = ion_histo_e_1st[0][traj_id]
        channel_occ_hist_x_shifted = [ion-1 for ion in channel_occ_hist_x[1:]]
        ax6.barh(channel_occ_hist_x_shifted, channel_occ_hist_y, align='center')

        # Plot 7 - L Coordination Timeseries Per Atom
        for ion_id in sorted(ion_ts_1st[0][traj_id].keys(), reverse=False):
            mpl.rcParams['axes.color_cycle'] = ['r', 'g', 'b', 'm', 'c']
            ion_ts_x = ion_ts_l_1st[1][traj_id][ion_id][::data_skip]
            ion_ts_y = ion_ts_l_1st[0][traj_id][ion_id][::data_skip]
            ion_ts_x_scale = [time_conv*pos for pos in ion_ts_x]
            ion_ts_y_flip = [1*pos for pos in ion_ts_y]
            ax7.plot(ion_ts_x_scale, ion_ts_y_flip)

        # Plot 8 - Coordination Histogram
        channel_occ_hist_x = ion_histo_l_1st[1][traj_id]
        channel_occ_hist_y = ion_histo_l_1st[0][traj_id]
        channel_occ_hist_x_shifted = [ion-1 for ion in channel_occ_hist_x[1:]]
        ax8.barh(channel_occ_hist_x_shifted, channel_occ_hist_y, align='center')

        # Plot 9 - E Coordination Timeseries Per Atom
        for ion_id in sorted(ion_ts_2nd[0][traj_id].keys(), reverse=False):
            mpl.rcParams['axes.color_cycle'] = ['r', 'g', 'b', 'm', 'c']
            ion_ts_x = ion_ts_e_2nd[1][traj_id][ion_id][::data_skip]
            ion_ts_y = ion_ts_e_2nd[0][traj_id][ion_id][::data_skip]
            ion_ts_x_scale = [time_conv*pos for pos in ion_ts_x]
            ion_ts_y_flip = [1*pos for pos in ion_ts_y]
            ax9.plot(ion_ts_x_scale, ion_ts_y_flip)

        # Plot 10 - Coordination Histogram
        channel_occ_hist_x = ion_histo_e_2nd[1][traj_id]
        channel_occ_hist_y = ion_histo_e_2nd[0][traj_id]
        channel_occ_hist_x_shifted = [ion-1 for ion in channel_occ_hist_x[1:]]
        ax10.barh(channel_occ_hist_x_shifted, channel_occ_hist_y, align='center')

        # Plot 11 - L Coordination Timeseries Per Atom
        for ion_id in sorted(ion_ts_2nd[0][traj_id].keys(), reverse=False):
            mpl.rcParams['axes.color_cycle'] = ['r', 'g', 'b', 'm', 'c']
            ion_ts_x = ion_ts_l_2nd[1][traj_id][ion_id][::data_skip]
            ion_ts_y = ion_ts_l_2nd[0][traj_id][ion_id][::data_skip]
            ion_ts_x_scale = [time_conv*pos for pos in ion_ts_x]
            ion_ts_y_flip = [1*pos for pos in ion_ts_y]
            ax11.plot(ion_ts_x_scale, ion_ts_y_flip)

        # Plot 12 - Coordination Histogram
        channel_occ_hist_x = ion_histo_l_2nd[1][traj_id]
        channel_occ_hist_y = ion_histo_l_2nd[0][traj_id]
        channel_occ_hist_x_shifted = [ion-1 for ion in channel_occ_hist_x[1:]]
        ax12.barh(channel_occ_hist_x_shifted, channel_occ_hist_y, align='center')

        # Plot 13 - E Coordination Timeseries Per Atom
        for ion_id in sorted(ion_ts_both[0][traj_id].keys(), reverse=False):
            mpl.rcParams['axes.color_cycle'] = ['r', 'g', 'b', 'm', 'c']
            ion_ts_x = ion_ts_e_both[1][traj_id][ion_id][::data_skip]
            ion_ts_y = ion_ts_e_both[0][traj_id][ion_id][::data_skip]
            ion_ts_x_scale = [time_conv*pos for pos in ion_ts_x]
            ion_ts_y_flip = [1*pos for pos in ion_ts_y]
            ax13.plot(ion_ts_x_scale, ion_ts_y_flip)

        # Plot 14 - Coordination Histogram
        channel_occ_hist_x = ion_histo_e_both[1][traj_id]
        channel_occ_hist_y = ion_histo_e_both[0][traj_id]
        channel_occ_hist_x_shifted = [ion-1 for ion in channel_occ_hist_x[1:]]
        ax14.barh(channel_occ_hist_x_shifted, channel_occ_hist_y, align='center')

        # Plot 15 - L Coordination Timeseries Per Atom
        for ion_id in sorted(ion_ts_both[0][traj_id].keys(), reverse=False):
            mpl.rcParams['axes.color_cycle'] = ['r', 'g', 'b', 'm', 'c']
            ion_ts_x = ion_ts_l_both[1][traj_id][ion_id][::data_skip]
            ion_ts_y = ion_ts_l_both[0][traj_id][ion_id][::data_skip]
            ion_ts_x_scale = [time_conv*pos for pos in ion_ts_x]
            ion_ts_y_flip = [1*pos for pos in ion_ts_y]
            ax15.plot(ion_ts_x_scale, ion_ts_y_flip)

        # Plot 16 - Coordination Histogram
        channel_occ_hist_x = ion_histo_l_both[1][traj_id]
        channel_occ_hist_y = ion_histo_l_both[0][traj_id]
        channel_occ_hist_x_shifted = [ion-1 for ion in channel_occ_hist_x[1:]]
        ax16.barh(channel_occ_hist_x_shifted, channel_occ_hist_y, align='center')

        # Plot 17 - Ion R Timeseries
        colors = ['r', 'g', 'b', 'm', 'c']
        for color_id, ion_id in enumerate(sorted(ion_ts_r[0][traj_id].keys(),
                                                 reverse=False)):
            ion_ts_x = ion_ts_r[1][traj_id][ion_id][::data_skip]
            ion_ts_y = ion_ts_r[0][traj_id][ion_id][::data_skip]
            ion_ts_x_scale = [time_conv*pos for pos in ion_ts_x]
            ion_ts_y_flip = [1*pos for pos in ion_ts_y]
            ax17.scatter(ion_ts_x_scale, ion_ts_y_flip, s=1, color=colors[color_id])

        # Plot 18 - Ion R Position Histogram
        for vals, ion_z in zip(ion_histo_r[0][traj_id], ion_histo_r[1][traj_id]):
            mpl.rcParams['axes.color_cycle'] = ['r', 'g', 'b', 'm', 'c']
            diff = (ion_z[1]-ion_z[0])/2.0
            ion_z_shifted = [1*(ion-diff) for ion in ion_z[1:]]
            vals_scaled = [val/20.0 for val in vals]
            ax18.plot(vals_scaled, ion_z_shifted)

        ax1.set_ylabel("Axial position (nm)")
        ax1.set_ylim([-0.7, 0.7])
        ax1.set_yticks([x/10.0 for x in range(-6,7,2)])
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
        ax5.set_ylabel("E\n1st")

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
        ax7.set_ylabel("L\n1st")

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
        ax9.set_ylabel("E\n2nd")

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
        ax11.set_ylabel("L\n2nd")
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
        ax13.set_ylabel("E\nBoth")

        ax14.grid(True)

        ax15.set_yticks(range(0,max_coord+1,2))
        ax15.set_ylim([-0.5,max_coord+0.5])
        ax15.set_ylabel("L\nBoth")
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
        ax17.set_ylim([0, 0.2])
        ax17.set_yticks([x/100.0 for x in range(0,20,10)])

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