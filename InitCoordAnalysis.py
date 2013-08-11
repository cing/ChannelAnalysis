#!/usr/bin/python

################################################################################
#
# Where it all begins. This script calls the entire data analysis pipeline
# based on a configuration input file.
#
# By Chris Ing, 2013 for Python 2.7
#
################################################################################
from ChannelAnalysis.CoordAnalysis import *
from ChannelAnalysis.RotamerAnalysis import *
from ChannelAnalysis.PoreAnalysis import *
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
        p = ConfigMap(cfgfile, 'plot_config')
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

        label_stats = label_statistics(data_f_states)

        '''
        print "Computing dunking populations"
        dunking_counts = rotamer_counter(data_rotamers,
                                         data_f_states,
                                         traj_col=int(r["traj_col"]))

        dunking_ts = compute_rotamer_vs_time(data_rotamers,
                                          data_f_states,
                                          traj_col=int(r["traj_col"]))
        #                                  prefix=r["out_prefix_rot"]+"rotamer")
        '''

        rotamer_2dhisto = compute_rotamer_rotamer_histogram(data_rotamers,
                                                           chi1_cols=chi1_cols,
                                                           chi2_cols=chi2_cols)

        plot_rotamer_2dhistograms(rotamer_2dhisto, state_dividers, label_stats,
                    prefix=p["out_prefix_plot"]+"figure8",
                    plot_title=p["fig_title_prefix"]+" N181 Chi1-Chi2 Histograms")

    raise " "

    # Extract SF residue time series information
    try:
        o = ConfigMap(cfgfile, 'oxygen_config')
        filenames_sf = ConfigMap(cfgfile, 'sf_residues').values()
        filenames_sf_pre = [o["input_file_prefix"] +
                             name for name in filenames_sf]

    except:
        print "No data for SF residues specified"
    else:
        print "Computing z deviation timeseries for sf chains"
        sf_processed = process_channelatoms(filenames_sf_pre,
                                         remove_frames=int(o["remove_frames"]))

        sf_ts = compute_atom_timeseries(sf_processed,
                                          int(o["sort_col"]),
                                          int(o["traj_col"]),
                                          col_skip=int(o["col_skip"]),
                                          num_cols=int(o["num_cols_per_ion"]),
                                          mean_shift=True)

    # Here we will attempt to extract as many oxygen positions in the
    # datafile as possible and append them to this list.
    allatom_histos = []

    # The same thing is true for oxygen atom Z timeseries data.
    alloxygen_ts = []

    # Extract residue 1 oxygen positions.
    try:
        o = ConfigMap(cfgfile, 'oxygen_config')
        filenames_res1 = ConfigMap(cfgfile, 'oxygen_residues1').values()
        filenames_res1_split = []
        for filename in filenames_res1:
            filenames_res1_split.extend(filename.split(','))
        filenames_res1_pre = [o["input_file_prefix"] +
                             name for name in filenames_res1_split]

    except:
        print "No oxygen data for residue 1 specified"
    else:
        print "Computing oxygen z histograms for residue 1 ring"
        oxygen_res1_processed = process_channelatoms(filenames_res1_pre,
                                         remove_frames=int(o["remove_frames"]))

        allatom_histos.append(compute_allatom_histogram(oxygen_res1_processed,
                                          int(o["sort_col"]),
                                          col_skip=int(o["col_skip"]),
                                          num_cols=int(o["num_cols_per_ion"])))

        alloxygen_ts.append(compute_atom_timeseries(oxygen_res1_processed,
                                          int(o["sort_col"]),
                                          int(o["traj_col"]),
                                          col_skip=int(o["col_skip"]),
                                          num_cols=int(o["num_cols_per_ion"])))

    # Extract residue 2 oxygen positions.
    try:
        o = ConfigMap(cfgfile, 'oxygen_config')
        filenames_res2 = ConfigMap(cfgfile, 'oxygen_residues2').values()
        filenames_res2_split = []
        for filename in filenames_res2:
            filenames_res2_split.extend(filename.split(','))
        filenames_res2_pre = [o["input_file_prefix"] +
                             name for name in filenames_res2_split]

    except:
        print "No oxygen data for residue 2 specified"
    else:
        print "Computing oxygen z histograms for residue 2 ring"
        oxygen_res2_processed = process_channelatoms(filenames_res2_pre,
                                         remove_frames=int(o["remove_frames"]))

        allatom_histos.append(compute_allatom_histogram(oxygen_res2_processed,
                                          int(o["sort_col"]),
                                          col_skip=int(o["col_skip"]),
                                          num_cols=int(o["num_cols_per_ion"])))

        alloxygen_ts.append(compute_atom_timeseries(oxygen_res2_processed,
                                          int(o["sort_col"]),
                                          int(o["traj_col"]),
                                          col_skip=int(o["col_skip"]),
                                          num_cols=int(o["num_cols_per_ion"])))

    # Extract residue 3 oxygen positions.
    try:
        o = ConfigMap(cfgfile, 'oxygen_config')
        filenames_res3 = ConfigMap(cfgfile, 'oxygen_residues3').values()
        filenames_res3_split = []
        for filename in filenames_res3:
            filenames_res3_split.extend(filename.split(','))
        filenames_res3_pre = [o["input_file_prefix"] +
                             name for name in filenames_res3_split]

    except:
        print "No oxygen data for residue 3 specified"
    else:
        print "Computing oxygen z histograms for residue 3 ring"
        oxygen_res3_processed = process_channelatoms(filenames_res3_pre,
                                         remove_frames=int(o["remove_frames"]))

        allatom_histos.append(compute_allatom_histogram(oxygen_res3_processed,
                                          int(o["sort_col"]),
                                          col_skip=int(o["col_skip"]),
                                          num_cols=int(o["num_cols_per_ion"])))


        alloxygen_ts.append(compute_atom_timeseries(oxygen_res3_processed,
                                          int(o["sort_col"]),
                                          int(o["traj_col"]),
                                          col_skip=int(o["col_skip"]),
                                          num_cols=int(o["num_cols_per_ion"])))

    # Extract residue 4 oxygen positions.
    try:
        o = ConfigMap(cfgfile, 'oxygen_config')
        filenames_res4 = ConfigMap(cfgfile, 'oxygen_residues4').values()
        filenames_res4_split = []
        for filename in filenames_res4:
            filenames_res4_split.extend(filename.split(','))
        filenames_res4_pre = [o["input_file_prefix"] +
                             name for name in filenames_res4_split]

    except:
        print "No oxygen data for residue 4 specified"
    else:
        print "Computing oxygen z histograms for residue 4 ring"
        oxygen_res4_processed = process_channelatoms(filenames_res4_pre,
                                         remove_frames=int(o["remove_frames"]))

        allatom_histos.append(compute_allatom_histogram(oxygen_res4_processed,
                                          int(o["sort_col"]),
                                          col_skip=int(o["col_skip"]),
                                          num_cols=int(o["num_cols_per_ion"])))

        alloxygen_ts.append(compute_atom_timeseries(oxygen_res4_processed,
                                          int(o["sort_col"]),
                                          int(o["traj_col"]),
                                          col_skip=int(o["col_skip"]),
                                          num_cols=int(o["num_cols_per_ion"])))

    # Extract residue 5 oxygen positions.
    try:
        o = ConfigMap(cfgfile, 'oxygen_config')
        filenames_res5 = ConfigMap(cfgfile, 'oxygen_residues5').values()
        filenames_res5_split = []
        for filename in filenames_res5:
            filenames_res5_split.extend(filename.split(','))
        filenames_res5_pre = [o["input_file_prefix"] +
                             name for name in filenames_res5_split]

    except:
        print "No oxygen data for residue 5 specified"
    else:
        print "Computing oxygen z histograms for residue 5 ring"
        oxygen_res5_processed = process_channelatoms(filenames_res5_pre,
                                         remove_frames=int(o["remove_frames"]))

        allatom_histos.append(compute_allatom_histogram(oxygen_res5_processed,
                                          int(o["sort_col"]),
                                          col_skip=int(o["col_skip"]),
                                          num_cols=int(o["num_cols_per_ion"])))

        alloxygen_ts.append(compute_atom_timeseries(oxygen_res5_processed,
                                          int(o["sort_col"]),
                                          int(o["traj_col"]),
                                          col_skip=int(o["col_skip"]),
                                          num_cols=int(o["num_cols_per_ion"])))

    print "Preparing plots for Oxygen Histograms"
    p = ConfigMap(cfgfile, 'plot_config')

    plot_oxy_timeseries(alloxygen_ts, sf_ts,
                     time_conv=float(p["time_conv"]),
                     prefix=p["out_prefix_plot"]+"figure7",
                     plot_title=p["fig_title_prefix"]+" E,L Oxygen Timeseries",
                     max_coord=int(p["max_coord"]),
                     max_length=int(p["max_length"]))

    # Extract all properties needed for 1st shell coordination studies.
    try:
        g = ConfigMap(cfgfile, 'coord_config')
        #print g
        filenames_1st = ConfigMap(cfgfile, 'coord_1st_input_files').values()
        filenames_1st_pre = [g["input_file_prefix"] +
                             name for name in filenames_1st]

        regex_pairs = sorted(ConfigMap(cfgfile, 'coord_state_label_regexs').iteritems())
        regexs = [value for key,value in regex_pairs]

        # This needs to be formatted from the comma-delimited string
        coord_cols = [int(col) for col in g["coordination_cols"].split(",")]
        sf_cols = [int(col) for col in g["selectivityf_cols"].split(",")]
        selectivityf_map = g["selectivityf_map"].split(",")
        selectivityf_code = "".join([res[0] for res in selectivityf_map])

        # These are plot variables that are unfortunately required at this time
        p = ConfigMap(cfgfile, 'plot_config')
        regexs_map = p["regex_map"].split(",")

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
                                 add_time=str2bool(g["add_time"]),
                                 padded=True)

        iondist_histo = compute_iondist_histograms(data_1st,
                                                int(g["sort_col"]),
                                       num_cols=int(g["num_cols_per_ion"]),
                                       occ_lower_cut=int(p["occ_lower_cut"]),
                                       occ_higher_cut=int(p["occ_higher_cut"]),
                                       histbins=150,
                                       prefix=g["out_prefix_1st"]+"1dhisto")


        plot_iondist_histograms(iondist_histo,
                                   prefix=p["out_prefix_plot"]+"figure6",
                plot_title=p["fig_title_prefix"]+" Pair Distance 1D Histograms")

        ionsplit_2dhisto = compute_ionsplit_2dhistograms(data_1st,
                                                int(g["sort_col"]),
                                       num_cols=int(g["num_cols_per_ion"]),
                                       occ_lower_cut=int(p["occ_lower_cut"]),
                                       occ_higher_cut=int(p["occ_higher_cut"]),
                                       histbins=150,
                                       prefix=g["out_prefix_1st"]+"2dhisto")

        plot_ionsplit_2dhistograms(ionsplit_2dhisto,
                                   prefix=p["out_prefix_plot"]+"figure3",
                             plot_title=p["fig_title_prefix"]+" 2D Histograms")

        print "Channel Occupancy and Figure 1 Calculations"
        counts_1st = occ_counter(data_1st,
                                 num_cols=int(g["num_cols_per_ion"]),
                                 traj_col=int(g["traj_col"]),
                                 sf_col=[])
        #                         prefix=g["out_prefix_1st"]+"chanocc")

        occ_1st = compute_occ_vs_time(data_1st,
                                      num_cols=int(g["num_cols_per_ion"]),
                                      traj_col=int(g["traj_col"]))
        #                              prefix=g["out_prefix_1st"]+"chanocc")

        ion_ts_1st = compute_ion_timeseries(data_1st,
                                           int(g["sort_col"]),
                                           int(g["traj_col"]),
                                           num_cols=int(g["num_cols_per_ion"]))

        ionsplit_histo = compute_ionsplit_histograms(data_1st,
                                                   int(g["sort_col"]),
                                          num_cols=int(g["num_cols_per_ion"]),
                                          prefix=g["out_prefix_1st"]+"1dhisto")

        allion_histo = compute_allion_histogram(data_1st,
                                                 int(g["sort_col"]),
                                          num_cols=int(g["num_cols_per_ion"]),
                                          prefix=g["out_prefix_1st"]+"1dhisto",
                                          normed=True)

        plot_ionsplit_histograms(ionsplit_histo, counts_1st,
                                 allion_histo, allatom_histos,
                                 prefix=p["out_prefix_plot"]+"figure2",
                                 plot_title=p["fig_title_prefix"]+" Ion Histograms")

        ion_histo_1st = compute_position_histograms(data_1st,
                                                    int(g["sort_col"]),
                                                    int(g["traj_col"]),
                                           num_cols=int(g["num_cols_per_ion"]))

        ion_ts_e_1st = compute_ion_timeseries(data_1st,
                                          int(g["sort_col"])+2,
                                          int(g["traj_col"]),
                                          num_cols=int(g["num_cols_per_ion"]))

        ion_histo_e_1st = compute_position_histograms(data_1st,
                                           int(g["sort_col"])+2,
                                           int(g["traj_col"]),
                                           num_cols=int(g["num_cols_per_ion"]),
                                           histmin=0, histmax=7,
                                           histbins=7)

        ion_ts_l_1st = compute_ion_timeseries(data_1st,
                                          int(g["sort_col"])+4,
                                          int(g["traj_col"]),
                                          num_cols=int(g["num_cols_per_ion"]))

        ion_histo_l_1st = compute_position_histograms(data_1st,
                                           int(g["sort_col"])+4,
                                           int(g["traj_col"]),
                                           num_cols=int(g["num_cols_per_ion"]),
                                           histmin=0, histmax=7,
                                           histbins=7)

        counts_sf_1st = occ_counter(data_1st,
                                    num_cols=int(g["num_cols_per_ion"]),
                                    traj_col=int(g["traj_col"]),
                                    sf_col=sf_cols)
        #                            prefix=g["out_prefix_1st"]+"sfocc")

        occ_sf_1st = compute_occ_vs_time(data_1st,
                                         num_cols=int(g["num_cols_per_ion"]),
                                         traj_col=int(g["traj_col"]),
                                         sf_col=sf_cols)
        #                                 prefix=g["out_prefix_1st"]+"sfocc")

        print "Writing multiple 1D histograms for ion coordination integers"
        coord_histo = compute_coord_histograms(data_1st,
                                          coord_cols=coord_cols,
                                          num_cols=int(g["num_cols_per_ion"]),
                                          sort_col=int(g["sort_col"]),
                                          prefix=g["out_prefix_1st"]+"1dhisto")

        print "Writing multiple 1D histograms for zero or non-zero coord"
        group_histo_1st = compute_group_coord_histograms(data_1st,
                                          sf_col=sf_cols,
                                          num_cols=int(g["num_cols_per_ion"]),
                                          sort_col=int(g["sort_col"]),
                                          prefix=g["out_prefix_1st"]+"1dhisto")

        print "Classifying sf columns as states with regular expressions"
        data_regex_1st = regex_columns(data_1st,
                                regex_strings=regexs,
                                num_cols=int(g["num_cols_per_ion"]),
                                sort_col=int(g["sort_col"]),
                                sort_cut=int(g["sort_cut"]),
                                sf_col=sf_cols,
                                max_ions=int(g["max_ions"]))

        print "Regex Macrostate Occupancy"
        counts_regex_1st = state_counter(data_1st, data_regex_1st,
                                         range(len(regexs)),
                                         traj_col=int(g["traj_col"]))

        print "Writing 1D histograms for regex coordination labels"
        histo_regex_1st = compute_regex_histograms(data_1st, data_regex_1st,
                                          traj_col=int(g["traj_col"]),
                                          num_cols=int(g["num_cols_per_ion"]),
                                          max_ions=int(g["max_ions"]),
                                          sort_col=int(g["sort_col"]),
                                          prefix=g["out_prefix_1st"]+"1dhisto")



        print "State Transition Counting"
        data_state_trans = state_transitions(data_1st, data_regex_1st,
                                             traj_col=int(g["traj_col"]))
        #print data_state_trans

        # Extract only the mean occupancies in percent (2nd last entry)
        data_f_mean_pop = counts_regex_1st[0]["MEAN"]
        print data_f_mean_pop

        print "State Transition Graph Building and Writing"
        #print data_state_trans

        state_draw_map = [(0.0, 0),(1.0, 0),(2.0, 0),(3.0, 0),
                              (0.5,-1),(1.5,-1),(2.5,-1),
                                   (1.0,-2),(2.0,-2),(3.0,-2),
                                       (1.5,-3),(2.5,-3),
                                            (2.0,-4)]

        #state_draw_map = [(0.0, 0),(1.0, 0),(0.5, -1)]

        data_state_trans_graph = build_macrostate_graph(data_state_trans,
                                                       state_draw_map,
                                                       pop_map=data_f_mean_pop)

        print "Macrostate Graph Writing"
        write_macrostate_graph(data_state_trans_graph, regexs_map,
                           plot_title=p["fig_title_prefix"] + " " +
                                      selectivityf_code + " " +
                                      "1st Shell Macrostates",
                           outfile=g["out_prefix_1st"] +
                                   "graph_regex_macro_1st_" +
                                   selectivityf_code + ".pdf")

        #print "Microstate Graph Writing"
        #write_microstate_graph(data_1st, data_regex_1st,
        #                  traj_col=int(g["traj_col"]),
        #                  outfile=g["out_prefix_1st"]+"graph_regex_micro.gexf")

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
                                 add_time=str2bool(g["add_time"]),
                                 padded=True)

        print "Channel Occupancy and Figure 1 Calculations"
        counts_2nd = occ_counter(data_2nd,
                                 num_cols=int(g["num_cols_per_ion"]),
                                 traj_col=int(g["traj_col"]),
                                 sf_col=[])
        #                         prefix=g["out_prefix_2nd"]+"chanocc")

        occ_2nd = compute_occ_vs_time(data_2nd,
                                      num_cols=int(g["num_cols_per_ion"]),
                                      traj_col=int(g["traj_col"]))
        #                              prefix=g["out_prefix_2nd"]+"chanocc")

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

        ion_histo_e_2nd = compute_position_histograms(data_2nd,
                                           int(g["sort_col"])+2,
                                           int(g["traj_col"]),
                                           num_cols=int(g["num_cols_per_ion"]),
                                           histmin=0, histmax=7,
                                           histbins=7)

        ion_ts_l_2nd = compute_ion_timeseries(data_2nd,
                                          int(g["sort_col"])+4,
                                          int(g["traj_col"]),
                                          num_cols=int(g["num_cols_per_ion"]))

        ion_histo_l_2nd = compute_position_histograms(data_2nd,
                                           int(g["sort_col"])+4,
                                           int(g["traj_col"]),
                                           num_cols=int(g["num_cols_per_ion"]),
                                           histmin=0, histmax=7,
                                           histbins=7)

        print "Writing multiple 1D histograms for zero or non-zero coord"
        group_histo_2nd = compute_group_coord_histograms(data_2nd,
                                          sf_col=sf_cols,
                                          num_cols=int(g["num_cols_per_ion"]),
                                          sort_col=int(g["sort_col"]),
                                          prefix=g["out_prefix_1st"]+"1dhisto")

    print "Made it passed 2nd"

    # Extract all properties needed for both shell coordination studies.
    try:
        g = ConfigMap(cfgfile, 'coord_config')
        #print g
        filenames_both = ConfigMap(cfgfile, 'coord_both_input_files').values()
        filenames_both_pre = [g["input_file_prefix"] +
                             name for name in filenames_both]

        regex_pairs = sorted(ConfigMap(cfgfile, 'coord_state_label_regexs').iteritems())
        regexs = [value for key,value in regex_pairs]

        # This needs to be formatted from the comma-delimited string
        coord_cols = [int(col) for col in g["coordination_cols"].split(",")]
        sf_cols = [int(col) for col in g["selectivityf_cols"].split(",")]
        selectivityf_map = g["selectivityf_map"].split(",")
        selectivityf_code = "".join([res[0] for res in selectivityf_map])

        # These are plot variables that are unfortunately required at this time
        p = ConfigMap(cfgfile, 'plot_config')
        regexs_map = p["regex_map"].split(",")

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
                                 add_time=str2bool(g["add_time"]),
                                 padded=True)

        print "Classifying sf columns as states with regular expressions"
        data_regex_both = regex_columns(data_both,
                                regex_strings=regexs,
                                num_cols=int(g["num_cols_per_ion"]),
                                sort_col=int(g["sort_col"]),
                                sort_cut=int(g["sort_cut"]),
                                sf_col=sf_cols,
                                max_ions=int(g["max_ions"]))

        print "Regex Macrostate Occupancy"
        counts_regex_both = state_counter(data_both, data_regex_both,
                                     range(len(regexs)),
                                     traj_col=int(g["traj_col"]))

        print "Writing 1D histograms for regex coordination labels"
        histo_regex_hoth = compute_regex_histograms(data_both, data_regex_both,
                                          traj_col=int(g["traj_col"]),
                                          num_cols=int(g["num_cols_per_ion"]),
                                          max_ions=int(g["max_ions"]),
                                          sort_col=int(g["sort_col"]),
                                          prefix=g["out_prefix_both"]+"1dhisto")

        print "State Transition Counting"
        data_state_trans = state_transitions(data_both, data_regex_both,
                                             traj_col=int(g["traj_col"]))

        # Extract only the mean occupancies in percent (2nd last entry)
        data_f_mean_pop = counts_regex_both[0]["MEAN"]
        #print data_f_mean_pop

        print "State Transition Graph Building and Writing"
        #print data_state_trans

        state_draw_map = [(0.0, 0),(1.0, 0),(2.0, 0),(3.0, 0),
                              (0.5,-1),(1.5,-1),(2.5,-1),
                                   (1.0,-2),(2.0,-2),(3.0,-2),
                                       (1.5,-3),(2.5,-3),
                                            (2.0,-4)]

        #state_draw_map = [(0.0, 0),(1.0, 0),(0.5, -1)]

        data_state_trans_graph = build_macrostate_graph(data_state_trans,
                                                       state_draw_map,
                                                       pop_map=data_f_mean_pop)

        print "Macrostate Graph Writing"
        write_macrostate_graph(data_state_trans_graph, regexs_map,
                           plot_title=p["fig_title_prefix"] + " " +
                                      selectivityf_code + " " +
                                      "Both Shell Macrostates",
                           outfile=g["out_prefix_both"] +
                                   "graph_regex_macro_both_" +
                                   selectivityf_code + ".pdf")

        #print "Microstate Graph Writing"
        #write_microstate_graph(data_1st, data_regex_1st,
        #                  traj_col=int(g["traj_col"]),
        #                  outfile=g["out_prefix_1st"]+"graph_regex_micro.gexf")

        print "State Transition w/ Intermediate Counting"
        print state_intermediate_transitions(data_both, data_regex_both,
                                             traj_col=int(g["traj_col"]))

        print "Channel Occupancy and Figure 1 Calculations"
        counts_both = occ_counter(data_both,
                                 num_cols=int(g["num_cols_per_ion"]),
                                 traj_col=int(g["traj_col"]),
                                 sf_col=[])
        #                         prefix=g["out_prefix_both"]+"chanocc")

        occ_both = compute_occ_vs_time(data_both,
                                      num_cols=int(g["num_cols_per_ion"]),
                                      traj_col=int(g["traj_col"]))
        #                              prefix=g["out_prefix_both"]+"chanocc")

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

        ion_histo_e_both = compute_position_histograms(data_both,
                                           int(g["sort_col"])+2,
                                           int(g["traj_col"]),
                                           num_cols=int(g["num_cols_per_ion"]),
                                           histmin=0, histmax=7,
                                           histbins=7)

        ion_ts_l_both = compute_ion_timeseries(data_both,
                                          int(g["sort_col"])+4,
                                          int(g["traj_col"]),
                                          num_cols=int(g["num_cols_per_ion"]))

        ion_histo_l_both = compute_position_histograms(data_both,
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

        print "Writing multiple 1D histograms for zero or non-zero coord"
        group_histo_both = compute_group_coord_histograms(data_both,
                                          sf_col=sf_cols,
                                          num_cols=int(g["num_cols_per_ion"]),
                                          sort_col=int(g["sort_col"]),
                                          prefix=g["out_prefix_1st"]+"1dhisto")

        print "Preparing plots for Stacked Histogram"

        plot_bindingmode_histograms(group_histo_1st,
                                    group_histo_2nd,
                                    group_histo_both,
                                 selectivityf_map=selectivityf_map,
                                 prefix=p["out_prefix_plot"]+"figure4",
                   plot_title=p["fig_title_prefix"]+" Binding Mode Histograms")

        plot_sf_w_2res_timeseries(occ_1st, counts_1st,
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
                           plot_title=p["fig_title_prefix"]+" Stacked Timeseries",
                           max_coord=int(p["max_coord"]),
                           max_length=int(p["max_length"]))


if __name__ == '__main__':
    cfgfile = ConfigParser()
    cfgfile.read(argv[1])
    main(cfgfile)
