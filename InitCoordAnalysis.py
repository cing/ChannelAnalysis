#!/usr/bin/python

################################################################################
#
# Where it all begins. This script calls the entire data analysis pipeline
# based on a configuration input file.
#
# By Chris Ing, 2013 for Python 2.7
#
################################################################################
from CoordAnalysis import *
from RotamerAnalysis import *
from CoordAnalysis import Coord_Plots
from sys import argv
from ConfigParser import ConfigParser

# A helper function to return values as dictionaries from the config file
def ConfigMap(cfgfile, section):
    dict1 = {}
    options = cfgfile.options(section)
    for option in options:
        try:
            dict1[option] = cfgfile.get(section, option)
            if dict1[option] == -1:
                DebugPrint("skip: %s" % option)
        except:
            print("exception on %s!" % option)
            dict1[option] = None
    return dict1

# A helper function to convert string expressions to boolean values.
def str2bool(v):
  return v.lower() in ("yes", "true", "t", "1")

# Extracts configuration data from the passed file and uses it to execute
# a number of calculations for plotting. In need of refactoring.
def main(cfgfile):

    # Extract all the properties needed for rotamer studies.
    try:
        r = ConfigMap(cfgfile, 'rotamer_config')
        #print r
        filenames_rot = ConfigMap(cfgfile, 'rotamer_input_files').values()
        filenames_rot_pre = [r["input_file_prefix"] +
                             name for name in filenames_rot]

        # This needs to be formatted from the comma-delimited string
        chi1_cols = [int(x1_val) for x1_val in r["x1_cols"].split(",")]
        chi2_cols = [int(x2_val) for x2_val in r["x2_cols"].split(",")]
        state_dividers = [int(div) for div in r["dividers"].split(",")]

    except:
        print "No rotamer data specified or no x1/x2 columns found"
    else:
        print "Performing rotamer data processing"
        # This is now dunking calculations.
        data_rotamers = process_rotamers(filenames=filenames_rot_pre,
                                         chi1_cols=chi1_cols,
                                         chi2_cols=chi2_cols,
                                         remove_frames=int(r["remove_frames"]),
                                         traj_col=int(r["traj_col"]))

        # Same thing but now we pass the SF column list.
        data_f_states = label_states(data_rotamers,
                                     chi2_cols,
                                     state_dividers)

        print "Computing dunking populations"
        dunking_counts = rotamer_counter(data_rotamers,
                                         data_f_states,
                                         traj_col=int(r["traj_col"]))

        dunking_ts = compute_rotamer_vs_time(data_rotamers,
                                          data_f_states,
                                          traj_col=int(r["traj_col"]),
                                          prefix=r["out_prefix_rot"]+"rotamer")

    # Extract all properties needed for 1st shell coordination studies.
    try:
        g = ConfigMap(cfgfile, 'coord_config')
        #print g
        filenames_1st = ConfigMap(cfgfile, 'coord_1st_input_files').values()
        filenames_1st_pre = [g["input_file_prefix"] +
                             name for name in filenames_1st]

        regexs = ConfigMap(cfgfile, 'coord_state_label_regexs').values()

        # This needs to be formatted from the comma-delimited string
        coord_cols = [int(col) for col in g["coordination_cols"].split(",")]
        sf_cols = [int(col) for col in g["selectivityf_cols"].split(",")]

    except:
        print "No 1st shell data specified or no coordination columns found"
    else:
        print "Sorting the ion coordination data for 1st shell"
        data_1st = process_input(filenames=filenames_1st_pre,
                                 num_cols=int(g["num_cols_per_ion"]),
                                 max_ions=int(g["max_ions"]),
                                 remove_frames=int(g["remove_frames"]),
                                 traj_col=int(g["traj_col"]),
                                 sort_col=int(g["sort_col"]),
                                 add_time=str2bool(g["add_time"]))

        print "Channel Occupancy and Figure 1 Calculations"
        counts_1st = occ_counter(data_1st,
                                 num_cols=int(g["num_cols_per_ion"]),
                                 traj_col=int(g["traj_col"]),
                                 sf_col=[],
                                 prefix=g["out_prefix_1st"]+"chanocc")

        occ_1st = compute_occ_vs_time(data_1st,
                                      num_cols=int(g["num_cols_per_ion"]),
                                      traj_col=int(g["traj_col"]),
                                      prefix=g["out_prefix_1st"]+"chanocc")

        ion_ts_1st = compute_ion_timeseries(data_1st,
                                           int(g["sort_col"]),
                                           int(g["traj_col"]),
                                           num_cols=int(g["num_cols_per_ion"]))

        ion_histo_1st = compute_position_histograms(data_1st,
                                                    int(g["sort_col"]),
                                                    int(g["traj_col"]),
                                           num_cols=int(g["num_cols_per_ion"]))

        ion_ts_e_1st = compute_ion_timeseries(data_1st,
                                          int(g["sort_col"])+2,
                                          int(g["traj_col"]),
                                          num_cols=int(g["num_cols_per_ion"]))

        ion_histo_e_1st = compute_position_histogram(data_1st,
                                           int(g["sort_col"])+2,
                                           int(g["traj_col"]),
                                           num_cols=int(g["num_cols_per_ion"]),
                                           histmin=0, histmax=7,
                                           histbins=7)

        ion_ts_l_1st = compute_ion_timeseries(data_1st,
                                          int(g["sort_col"])+4,
                                          int(g["traj_col"]),
                                          num_cols=int(g["num_cols_per_ion"]))

        ion_histo_l_1st = compute_position_histogram(data_1st,
                                           int(g["sort_col"])+4,
                                           int(g["traj_col"]),
                                           num_cols=int(g["num_cols_per_ion"]),
                                           histmin=0, histmax=7,
                                           histbins=7)

        counts_sf_1st = occ_counter(data_1st,
                                    num_cols=int(g["num_cols_per_ion"]),
                                    traj_col=int(g["traj_col"]),
                                    sf_col=sf_cols,
                                    prefix=g["out_prefix_1st"]+"sfocc")

        occ_sf_1st = compute_occ_vs_time(data_1st,
                                         num_cols=int(g["num_cols_per_ion"]),
                                         traj_col=int(g["traj_col"]),
                                         sf_col=sf_cols,
                                         prefix=g["out_prefix_1st"]+"sfocc")

        print "Writing multiple 1D histograms for ion coordination integers"
        coord_histo = compute_coord_histograms(data_1st,
                                          coord_cols=coord_cols,
                                          num_cols=int(g["num_cols_per_ion"]),
                                          sort_col=int(g["sort_col"]),
                                          prefix=g["out_prefix_1st"]+"1dhisto")

        print "Writing multiple 1D histograms for zero or non-zero coord"
        group_histo = compute_group_coord_histograms(data_1st,
                                          sf_col=sf_cols,
                                          num_cols=int(g["num_cols_per_ion"]),
                                          sort_col=int(g["sort_col"]),
                                          prefix=g["out_prefix_1st"]+"1dhisto")

        print "Padding the ion coordination data for 1st shell"
        data_1st = process_input(filenames=filenames_1st_pre,
                                 num_cols=int(g["num_cols_per_ion"]),
                                 max_ions=int(g["max_ions"]),
                                 remove_frames=int(g["remove_frames"]),
                                 traj_col=int(g["traj_col"]),
                                 sort_col=int(g["sort_col"]),
                                 add_time=str2bool(g["add_time"]),
                                 padded=True)

        print "Classifying sf columns as states with regular expressions"
        data_regex_1st = regex_columns(data_1st,
                                regex_strings=regexs,
                                num_cols=int(g["num_cols_per_ion"]),
                                sort_col=int(g["sort_col"]),
                                sort_cut=int(g["sort_cut"]),
                                sf_col=sf_cols,
                                max_ions=int(g["max_ions"]))

        print "Regex Macrostate Occupancy"
        counts_regex = regex_counter(data_1st, data_regex_1st,
                                     traj_col=int(g["traj_col"]))

        print "Writing 1D histograms for regex coordination labels"
        histo_regex = compute_regex_histograms(data_1st, data_regex_1st,
                                          traj_col=int(g["traj_col"]),
                                          num_cols=int(g["num_cols_per_ion"]),
                                          max_ions=int(g["max_ions"]),
                                          sort_col=int(g["sort_col"]),
                                          prefix=g["out_prefix_1st"]+"1dhisto")

        print "State Transition Counting"
        data_state_trans = state_transitions(data_1st, data_regex_1st,
                                             traj_col=int(g["traj_col"]))
        print data_state_trans

        # Extract only the mean occupancies in percent (2nd last entry)
        data_f_mean_pop = counts_regex[0]["MEAN"]
        print data_f_mean_pop

        print "State Transition Graph Building and Writing"
        #print data_state_trans
        '''
        state_draw_map = [(0.0, 0),(1.0, 0),(2.0, 0),(3.0, 0),
                              (0.5,-1),(1.5,-1),(2.5,-1),
                                   (1.0,-2),(2.0,-2),(3.0,-2),
                                       (1.5,-3),(2.5,-3),
                                            (2.0,-4)]
        '''
        state_draw_map = [(0.0, 0),(1.0, 0),(0.5, -1)]

        data_state_trans_graph = build_macrostate_graph(data_state_trans,
                                                       state_draw_map,
                                                       pop_map=data_f_mean_pop)

        print "Macrostate Graph Writing"
        write_macrostate_graph(data_state_trans_graph,
                           outfile=g["out_prefix_1st"]+"graph_regex_macro.pdf")

        print "Microstate Graph Writing"
        write_microstate_graph(data_1st, data_regex_1st,
                          traj_col=int(g["traj_col"]),
                          outfile=g["out_prefix_1st"]+"graph_regex_micro.gexf")

        print "State Transition w/ Intermediate Counting"
        print state_intermediate_transitions(data_1st, data_regex_1st,
                                             traj_col=int(g["traj_col"]))


    # Extract all properties needed for 2nd shell coordination studies.
    try:
        g = ConfigMap(cfgfile, 'coord_config')
        #print g
        filenames_2nd = ConfigMap(cfgfile, 'coord_2nd_input_files').values()
        filenames_2nd_pre = [g["input_file_prefix"] +
                             name for name in filenames_2nd]
    except:
        print "No 2nd shell data specified"
    else:
        print "Sorting the ion coordination data for 2nd shell"
        data_2nd = process_input(filenames=filenames_2nd_pre,
                                 num_cols=int(g["num_cols_per_ion"]),
                                 max_ions=int(g["max_ions"]),
                                 remove_frames=int(g["remove_frames"]),
                                 traj_col=int(g["traj_col"]),
                                 sort_col=int(g["sort_col"]),
                                 add_time=str2bool(g["add_time"]))

        print "Channel Occupancy and Figure 1 Calculations"
        counts_2nd = occ_counter(data_2nd,
                                 num_cols=int(g["num_cols_per_ion"]),
                                 traj_col=int(g["traj_col"]),
                                 sf_col=[],
                                 prefix=g["out_prefix_2nd"]+"chanocc")

        occ_2nd = compute_occ_vs_time(data_2nd,
                                      num_cols=int(g["num_cols_per_ion"]),
                                      traj_col=int(g["traj_col"]),
                                      prefix=g["out_prefix_2nd"]+"chanocc")

        ion_ts_2nd = compute_ion_timeseries(data_2nd,
                                           int(g["sort_col"]),
                                           int(g["traj_col"]),
                                           num_cols=int(g["num_cols_per_ion"]))

        ion_histo_2nd = compute_position_histograms(data_2nd,
                                                    int(g["sort_col"]),
                                                    int(g["traj_col"]),
                                           num_cols=int(g["num_cols_per_ion"]))

        ion_ts_e_2nd = compute_ion_timeseries(data_2nd,
                                          int(g["sort_col"])+2,
                                          int(g["traj_col"]),
                                          num_cols=int(g["num_cols_per_ion"]))

        ion_histo_e_2nd = compute_position_histogram(data_2nd,
                                           int(g["sort_col"])+2,
                                           int(g["traj_col"]),
                                           num_cols=int(g["num_cols_per_ion"]),
                                           histmin=0, histmax=7,
                                           histbins=7)

        ion_ts_l_2nd = compute_ion_timeseries(data_2nd,
                                          int(g["sort_col"])+4,
                                          int(g["traj_col"]),
                                          num_cols=int(g["num_cols_per_ion"]))

        ion_histo_l_2nd = compute_position_histogram(data_2nd,
                                           int(g["sort_col"])+4,
                                           int(g["traj_col"]),
                                           num_cols=int(g["num_cols_per_ion"]),
                                           histmin=0, histmax=7,
                                           histbins=7)

    # Extract all properties needed for both shell coordination studies.
    try:
        g = ConfigMap(cfgfile, 'coord_config')
        #print g
        filenames_both = ConfigMap(cfgfile, 'coord_both_input_files').values()
        filenames_both_pre = [g["input_file_prefix"] +
                             name for name in filenames_both]
    except:
        print "No both shell data specified"
    else:
        print "Sorting the ion coordination data for both shells"
        data_both = process_input(filenames=filenames_both_pre,
                                 num_cols=int(g["num_cols_per_ion"]),
                                 max_ions=int(g["max_ions"]),
                                 remove_frames=int(g["remove_frames"]),
                                 traj_col=int(g["traj_col"]),
                                 sort_col=int(g["sort_col"]),
                                 add_time=str2bool(g["add_time"]))

        print "Channel Occupancy and Figure 1 Calculations"
        counts_both = occ_counter(data_both,
                                 num_cols=int(g["num_cols_per_ion"]),
                                 traj_col=int(g["traj_col"]),
                                 sf_col=[],
                                 prefix=g["out_prefix_both"]+"chanocc")

        occ_both = compute_occ_vs_time(data_both,
                                      num_cols=int(g["num_cols_per_ion"]),
                                      traj_col=int(g["traj_col"]),
                                      prefix=g["out_prefix_both"]+"chanocc")

        ion_ts_both = compute_ion_timeseries(data_both,
                                           int(g["sort_col"]),
                                           int(g["traj_col"]),
                                           num_cols=int(g["num_cols_per_ion"]))

        ion_histo_both = compute_position_histograms(data_both,
                                                    int(g["sort_col"]),
                                                    int(g["traj_col"]),
                                           num_cols=int(g["num_cols_per_ion"]))

        ion_ts_e_both = compute_ion_timeseries(data_both,
                                          int(g["sort_col"])+2,
                                          int(g["traj_col"]),
                                          num_cols=int(g["num_cols_per_ion"]))

        ion_histo_e_both = compute_position_histogram(data_both,
                                           int(g["sort_col"])+2,
                                           int(g["traj_col"]),
                                           num_cols=int(g["num_cols_per_ion"]),
                                           histmin=0, histmax=7,
                                           histbins=7)

        ion_ts_l_both = compute_ion_timeseries(data_both,
                                          int(g["sort_col"])+4,
                                          int(g["traj_col"]),
                                          num_cols=int(g["num_cols_per_ion"]))

        ion_histo_l_both = compute_position_histogram(data_both,
                                           int(g["sort_col"])+4,
                                           int(g["traj_col"]),
                                           num_cols=int(g["num_cols_per_ion"]),
                                           histmin=0, histmax=7,
                                           histbins=7)

        ion_ts_r = compute_ion_sqrt_timeseries(data_both,
                                          [int(g["sort_col"])-2,
                                           int(g["sort_col"])-1],
                                          int(g["traj_col"]),
                                          num_cols=int(g["num_cols_per_ion"]))

        ion_histo_r = compute_position_sqrt_histograms(data_both,
                                           [int(g["sort_col"])-2,
                                           int(g["sort_col"])-1],
                                           int(g["traj_col"]),
                                           num_cols=int(g["num_cols_per_ion"]))

        '''
        print "Preparing plots for Stacked Histogram"
        p = ConfigMap(cfgfile, 'plot_config')
        plot_stacked_timeseries(occ_1st, counts_1st,
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
                                time_conv=float(p["time_conv"]),
                                prefix=p["out_prefix_plot"]+"figure1",
                                plot_title=p["fig1_title"],
                                max_coord=int(p["max_coord"]),
                                max_length=int(p["max_length"]))
        '''


if __name__ == '__main__':
    cfgfile = ConfigParser()
    cfgfile.read(argv[1])
    main(cfgfile)
