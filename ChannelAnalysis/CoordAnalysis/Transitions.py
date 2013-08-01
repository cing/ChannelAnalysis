#!/usr/bin/python

################################################################################
#
# This script parses a state stream and detects transitions that
# satisfy a particular regular expression for the state in start,
# N intermediates, and an end state. The script can exclude short-lived
# states and function with multiple trajectory datafiles as input.
# This script builds on the principles of the
# Dwell_Counter_State_Labels_VariableMem.py and Dwell_Counter_State_Labels.py
# scripts.
#
# Example: Given the input state stream where we intend to detect 1 2 3:
#          3
#          1
#          1
#          2
#          2
#          2
#          4
#          2
#          3
#          3
#          3
#          2
#          This script would still detect a 1 2 3 transition if we had
#          enforced a label_memory of 3 and a threshhold >= 0.8
#
# By Chris Ing, 2013 for Python 2.7
#
################################################################################
from ChannelAnalysis.CoordAnalysis.Preprocessor import *
from ChannelAnalysis.CoordAnalysis.Grouping import *
from re import match
from argparse import ArgumentParser
from collections import defaultdict
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

# This function will extract all transitions (including fast recrossings)
# and return a list with key of type "regex_label1"-"regex_label2"
# and the value is the number of those transitions observed.
def state_transitions(data_floats, data_regex,
                      time_col=0, time_increment=1,
                      traj_col=11, verbose=False):

    transition_totals = defaultdict(int)

    # the data_regex datatype is a list of tuples: (state_label, regex_int)
    data_regex_ids = [data[1] for data in data_regex]

    # Extract the useful data out of these large lists
    data_paired = []
    for line, regid in zip(data_floats,data_regex_ids):
        traj = line[traj_col]
        time = line[time_col]
        data_paired.append([time,traj,regid])

    # n[0] is time, n[1] is trajectory, n[2] is regex_state
    for n, nplus1 in window(data_paired,2):
        if verbose:
            print (n, nplus1, n[1] == nplus1[1],
                  nplus1[0] - n[0] == time_increment)
        # Verify there is no missing timestep and that the run number
        # has not changed, also verify that it's a transition.
        if (n[1] == nplus1[1] and
            nplus1[0] - n[0] == time_increment and
            n[2] != nplus1[2]):
            #print N[1], Nplus1[1]
            transition_totals[str(n[2])+"-"+str(nplus1[2])] += 1

    return transition_totals.items()

# This function writes out a gexf graph file of the entire microstate network
# which can then be visualized in Gephi.
def write_microstate_graph(data_floats, data_regex, time_col=0,
                           time_increment=1, traj_col=11,
                           outfile="graph_regex_micro.gexf"):
    try:
        import networkx as nx
    except ImportError:
        print "This function requires the networkx package for graph building"
    else:
        micro_graph = nx.DiGraph()

        # the data_regex datatype is a list of tuples: (state_label, regex_int)
        state_labels = [data[0] for data in data_regex]
        regex_ids = [data[1] for data in data_regex]

        # Extract the useful data out of these large lists
        data_paired = []
        for line, regid, state in zip(data_floats, regex_ids, state_labels):
            traj = line[traj_col]
            time = line[time_col]
            data_paired.append([time,traj,regid,state])

        # Build nodes with dwell times data
        for label in state_labels:
            if micro_graph.has_node(label):
                micro_graph.node[label]["time"] += 1
            else:
                micro_graph.add_node(label, time=1)

        # Now loop over pairs to add connectivity.
        # n[0] is time, n[1] is trajectory, n[2] is regex_state, n[3] is label
        for n, nplus1 in window(data_paired,2):
            # Verify there is no missing timestep and that the run number
            # has not changed, also verify that it's a transition.
            if (n[1] == nplus1[1] and
                nplus1[0] - n[0] == time_increment and
                n[3] != nplus1[3]):
                if micro_graph.has_edge(n[3], nplus1[3]):
                    micro_graph.edge[n[3]][nplus1[3]]["weight"] += 1
                else:
                    micro_graph.add_edge(n[3], nplus1[3], weight=1)

        nx.write_gexf(micro_graph, outfile)

    return micro_graph

# This function processes a state stream and interally turns it into a
# list of the format: time traj state_id time_spent_in_that_state
#
# Using this format it is trivial to exclude low occupancy intermediates
# and extract N-state permeation events. It is up to the user to interpret
# which of the N-state events are biologically revelant and indicate a
# permeation event. The intermediates parameter adjusts the N variable
# although typically 1 is sufficient.
# The datatype returned is a list of tuples where the first tuple index
# is the transition id given by a sequence of hypen separated regex ids
# and the value is the transition count.
def state_intermediate_transitions(data_floats, data_regex,
                      time_col=0, time_increment=1, intermediates=1,
                      traj_col=11, cutoff=1, verbose=False):

    # the data_regex datatype is a list of tuples: (state_label, regex_int)
    data_regex_ids = [data[1] for data in data_regex]

    # Extract the useful data out of these large lists with the traj
    # as the key and the list as the timesteps.
    data_paired_by_traj = defaultdict(list)
    for line, regid in zip(data_floats,data_regex_ids):
        traj = line[traj_col]
        time = line[time_col]
        data_paired_by_traj[traj].append([time,traj,regid])

    data_grouped = []
    # This time we'll group the state stream into a list of lists grouped by
    # regex_id with the magic of itertools.groupby:
    # http://docs.python.org/2/library/itertools.html#itertools.groupby
    for traj, timesteps in data_paired_by_traj.iteritems():
        for key, group in groupby(timesteps, key=lambda col: col[2]):
            group_list = list(group)
            group_len = len(group_list)*time_increment

            # Here we can impose a cutoff that removes low-population
            # states in between large intermediates.
            if group_len > cutoff:
                data_grouped.append(group_list[0]+[group_len])
            #print group_list[0]+[group_len]

    # Any cutoff value other than 0 will create gaps in the stream
    # and these need to be corrected. The following block of code
    # will look for gaps and adjust state lifetimes accordingly.
    data_collapsed = []
    # Make sure we have at least 1 group and then extract it's regex ID.
    assert len(data_grouped) > 1
    previous_group = data_grouped[0]
    start_group = data_grouped[0]
    temp_dwell = 0
    for current_group in data_grouped:
        # If there is a change of traj num since the past step
        # then add the final dwell_time value before resetting
        if previous_group[1] != current_group[1]:
            temp_dwell += previous_group[3]
            # We only take the first 3 columns because we're rewriting
            # the last column.
            data_collapsed.append(start_group[:3]+[temp_dwell])
            temp_dwell = 0
            start_group = current_group
        # If there is a change of id since the past step
        # use the difference in time values to correct the dwell time
        elif previous_group[2] != current_group[2]:
            temp_dwell += current_group[0]-previous_group[0]
            data_collapsed.append(start_group[:3]+[temp_dwell])
            temp_dwell = 0
            start_group = current_group
        else:
            temp_dwell += current_group[0]-previous_group[0]
            #print current_group[0], temp_dwell
        previous_group = current_group

    temp_dwell += previous_group[0]-start_group[0]
    data_collapsed.append(start_group[:3]+[temp_dwell])

    # This dictionary is returned by the function and summarizes the observed
    # transitions.
    transition_totals = defaultdict(int)
    # intermediates+2 is used because we have start and end states to consider
    for group_window in window(data_collapsed,intermediates+2):
        # n[0] is time, n[1] is traj_id, n[2] is regex_id, n[3] is dwell_time
        group_regex_ids = [str(temp_group[2]) for temp_group in group_window]
        # Confirm that a real transition ocurred oppossed to a back-crossing
        if len(group_regex_ids) == len(set(group_regex_ids)):
            transition_totals["-".join(group_regex_ids)] += 1

    return transition_totals.items()

# This function will use networkx to build a graph object of nodes, edges.
# The state_list is a list of tuples where each tuple is a transition ID
# and a transition count. The transition ID is a "delimited" list of state_ids
# starting from zero. All possible states in the state_list
# array need to be placed inside a diamond using the grid_map. grid_map is a
# list of strings that represent positions on a 2D lattice. Anything in
# the 2D lattice without a corresponding grid_map entry will be deleted.
def build_macrostate_graph(state_edges, grid_map, pop_map=None, delim="-",
                           count_cut=10):
    try:
        import networkx as nx
    except ImportError:
        print "This function requires the networkx package for graph building"
    else:
        # Identify the ID for unclassified data
        skip_id = max([int(a) for i,c in state_edges for a in i.split(delim)])

        macro_graph = nx.DiGraph()
        for trans_id, count in state_edges:
            nodes = trans_id.split(delim)
            for node in nodes:
                # exclude the skip_id computed above
                if int(node) != skip_id:
                    assert int(node) < len(grid_map)
                    if not macro_graph.has_node(node):
                        # pop_map stores the population of each regex state
                        if pop_map is not None:
                            trunc_pop = "%.1f" % (100.*pop_map[int(node)])
                            macro_graph.add_node(node,
                                                   pos=grid_map[int(node)],
                                                   weight=trunc_pop)
                        else:
                            macro_graph.add_node(node,
                                                   pos=grid_map[int(node)])

            if count > count_cut:
                if (int(nodes[0]) != skip_id) and (int(nodes[-1]) != skip_id):
                    macro_graph.add_edge(nodes[0],nodes[-1], weight=count)

        return macro_graph

def write_macrostate_graph(macro_graph, regexs_map=None,
                           plot_title="Macrostate Graph",
                           outfile="graph_regex_macro.pdf"):
    try:
        import networkx as nx
        import matplotlib.pyplot as plt
    except ImportError:
        print "This function requires networkx/matplotlib for graph plotting"
    else:
        # Initialize a new plot
        fig = plt.figure()

        pos=nx.get_node_attributes(macro_graph,'pos')
        elabels=nx.get_edge_attributes(macro_graph,'weight')
        nlabels=nx.get_node_attributes(macro_graph,'weight')

        # Here we're going to add the state populations to the node labels
        if len(nlabels) > 0:
            for label, pop in nlabels.iteritems():
                if regexs_map != None:
                    nlabels[label] = regexs_map[int(label)] + \
                                     "\n(" + nlabels[label] + ")"
                else:
                    nlabels[label] = label + "\n(" + nlabels[label] + ")"

        # Since networkx graph drawing actually kinda stinks for DiGraph's
        # I have to make joint labels for forward and reverse edges.
        uniq_labels={}
        for label, count in elabels.iteritems():
            rev_label = label[::-1]
            #print label, rev_label, count
            if (label not in uniq_labels) and (rev_label not in uniq_labels):
                uniq_labels[label] = str(count)
            elif (label not in uniq_labels) and (rev_label in uniq_labels):
                uniq_labels[rev_label] += "/" + str(count)
            else:
                uniq_labels[label] += "/" + str(count)

        nx.draw_networkx(macro_graph, pos,
                         with_labels=True,
                         node_shape='s',
                         node_size=600,
                         labels=nlabels,
                         font_size=5)
        nx.draw_networkx_edge_labels(macro_graph, pos,
                                     edge_labels=uniq_labels,
                                     font_size=5,
                                     label_pos=0.40)
        plt.axis('off')
        plt.subplots_adjust(hspace = 0.1, wspace = 0.02,
                            left = 0.1, bottom = 0.08,
                            right = 0.95, top = 0.75)
        plt.suptitle(plot_title)
        plt.savefig(outfile)

if __name__ == '__main__':
    parser = ArgumentParser(
    description='This script takes regex expressions for a state label\
    and outputs how many of your states are classified by that label\
    as well as the population of those states in the dataset')

    parser.add_argument(
    '-f', dest='filenames', type=str, nargs="+", required=True,
    help='a filename of coordination data from MDAnalysis trajectory data')
    parser.add_argument(
    '-m', dest='max_ions', type=int, required=True,
    help='the maximum number of ions in the channel to consider')
    parser.add_argument(
    '-c', dest='num_cols', type=int, default=13,
    help='the number of columns per ion in the input')
    parser.add_argument(
    '-remove', dest='remove_frames', type=int, default=0,
    help='this is a number of frames to remove from the start of the data')
    parser.add_argument(
    '-s', dest='sort_col', type=int, default=3,
    help='a zero inclusive column number to sort your row on, typically x,y,z')
    parser.add_argument(
    '-sc', dest='sort_cut', type=float, default=0.0,
    help='a value on the sort_col range to classify zero coordinated data')
    parser.add_argument(
    '-sf', dest='sf_col', type=int, nargs="+", default=[5,6],
    help='the coordination integer columns that define the selectivity filter')
    parser.add_argument(
    '-n', dest='num_ions_map', type=int, nargs="+", default=[1,2,2,3,0],
    help='list of integer ion counts in SF for each regex value + 1 for extra')
    parser.add_argument(
    '-t', dest='traj_col', type=int, default=11,
    help='a zero inclusive column number that contains the run number')
    parser.add_argument(
    '-o', dest='outfile', type=str, default=None,
    help='the file to output the sorted padding output of all input files')
    parser.add_argument(
    '--addtime', dest='add_time', action="store_true", default=False,
    help='an optional argument to add time columns to each ion grouping')
    parser.add_argument(
    '-i', dest='regex', type=str, nargs="+", required=True,
    help='a list of regex values in quotes')
    parser.add_argument(
    '-timecol', dest='time_col', type=int, default=0,
    help='a zero inclusive column number with the frame/timestep number')
    parser.add_argument(
    '-dt', dest='time_increment', type=int, default=1,
    help='the difference in the time column between steps')
    args = parser.parse_args()

    data_f_padded = process_input(filenames=args.filenames,
                                          num_cols=args.num_cols,
                                          max_ions=args.max_ions,
                                          remove_frames=args.remove_frames,
                                          traj_col=args.traj_col,
                                          sort_col=args.sort_col,
                                          add_time=args.add_time,
                                          padded=True)

    data_f_regex = regex_columns(data_f_padded, regex_strings=args.regex,
                                num_cols=args.num_cols,
                                sort_col=args.sort_col,
                                sort_cut=args.sort_cut,
                                sf_col=args.sf_col,
                                max_ions=args.max_ions)

    print "State Transition Counting"
    data_state_trans = state_transitions(data_f_padded, data_f_regex,
                            time_col=args.time_col,
                            time_increment=args.time_increment,
                            traj_col=args.traj_col)

    print "Computing Regex State Occupancies for macrostate graph"
    data_f_occupancy = regex_counter(data_f_padded, data_f_regex,
                        num_ions_map=args.num_ions_map,
                        traj_col=args.traj_col)
    # We want to extract only the mean occupancies in percent (2nd last entry)
    data_f_mean_pop = data_f_occupancy[-2][1:]

    print "State Transition Graph Building"
    #print data_state_trans
    state_draw_map = [(0.0, 0),(1.0, 0),(2.0, 0),(3.0, 0),
                          (0.5,-1),(1.5,-1),(2.5,-1),
                               (1.0,-2),(2.0,-2),(3.0,-2),
                                   (1.5,-3),(2.5,-3),
                                        (2.0,-4)]

    data_state_trans_graph = build_macrostate_graph(data_state_trans,
                                                    state_draw_map,
                                                    pop_map=data_f_mean_pop)

    print "Macrostate Graph Writing"
    write_macrostate_graph(data_state_trans_graph)

    print "Microstate Graph Writing"
    write_microstate_graph(data_f_padded, data_f_regex,
                           time_col=args.time_col,
                           time_increment=args.time_increment,
                           traj_col=args.traj_col)

    '''
    print "State Transition w/ Intermediate Counting"
    print state_intermediate_transitions(data_f_padded, data_f_regex,
                                         time_col=args.time_col,
                                         time_increment=args.time_increment,
                                         traj_col=args.traj_col)
    '''
