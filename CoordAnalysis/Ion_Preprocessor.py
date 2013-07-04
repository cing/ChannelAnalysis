#!/usr/bin/python
###############################################################################
#
# This script prepares multiple files of coordination output from MDAnalysis
# coordination output scripts in order to facilitate numerous 1D distributions,
# graph building, and state splitting. This script performs a variety of
# simple processing functions including adding columns, padding the data
# with 0's and sorting the columns by a particular column (namely z value).
#
# Example: For 13/26/39/...-column data with this type (describe elsewhere):
#          1.0 -0.13 -0.193 0.522 0.0 0.0 0.0 0.0 0.0 2.0 9.0 2.0 1748.0
#          2.0 -0.124 -0.013 0.662 0.0 0.0 1.0 0.0 0.0 2.0 8.0 2.0 1748.0
#
#          The following command will load two datafiles, normalize
#          the column length to 39 and remove the first 2000 lines:
#
#          python Ion_Preprocessor.py -f f1.out f2.out -m 3 -c 13 -remove 2000
#
#
# This script typically is the second step in a larger analysis pipeline.
# As far as script history, it's a combination of Row_Add_Time.py,
# Row_Pad_Columns.py, Row_Sorted_OnZKeepAllTags.py, ion_count_splitter.py.
#
# By Chris Ing, 2013 for Python 2.7
#
###############################################################################
from argparse import ArgumentParser
from re import match

#a great helper function to iterate over chunks of a list
def chunker(seq, size):
    return (seq[pos:pos + size] for pos in xrange(0, len(seq), size))

# This script takes a list of lists containing floats
# and adds time. Some datasets do not require this step.
def add_time_column(data_floats, num_cols=13, verbose=False):

    data_output = []
    for line in data_floats:
        non_time_columns = line[1:]

        # We chunk the non-time columns into ion groupings that
        # were detected at that timestep (-1 because there is no time col)
        temp_line = []
        for ion in chunker(non_time_columns,num_cols-1):
            temp_line.append(line[0])
            temp_line.extend(ion)
            if verbose:
                print line[0], " ".join([str(ion_col) for ion_col in ion]),

        # Here we're building a series of temp lines and appending
        # one for each timestep.
        data_output.append(temp_line)
        if verbose:
            print "\n",

    return data_output

# This is a script that produces a regex classification from the input data
# instead of recomputing it multiple times throughout the script for things
# like transition counting. It requires sorted or sorted/padded data
# that are returned from sort_columns and pad_columns functions.
def regex_columns(data_floats, regex_strings, num_cols=13,
                  pad_col=4, sf_col=[5,6], sort_col=3, sort_cut=0.0,
                  max_ions=3, max_coord=9, prefix=None):

    data_output = []

    # File streams for splitting purposes if the prefix flag is set
    if prefix != None:
        count_files={}
        for regex_id in range(len(regex_strings)+1):
            count_files[regex_id] = open(prefix+"_split"+str(regex_id),"w")

    for line in data_floats:
        temp_label = []
        num_ions = 0
        for ion in chunker(line,num_cols):

            # This functionality was implemented upon discovering that
            # zero coordinated ions can exist in either the central cavity
            # or the extracellular region of a channel. For this case,
            # a sort_value cutoff must be used to determine where the ion
            # actually is! Note, look at your SF bound ion distributions
            # to see if sort_cut is set correctly.
            sort_val = ion[sort_col]
            all_zeros = all([ion[col]==0.0 for col in sf_col])

            # I need to code a bit of logic for pre-processed data.
            if ion[pad_col] != "-":
                if (sort_val > sort_cut) and all_zeros:
                    temp_label.extend("+"*len(sf_col))
                else:
                    temp_label.extend([int(ion[col]) for col in sf_col])
            else:
                temp_label.extend("-"*len(sf_col))
            num_ions += 1

        # Here's a fix for when coordination integer counts are too large
        # and it ruins the fixed number of digit state label paradigm that
        # is critical for regex matching.
        for digit_index, digit in enumerate(temp_label):
            if digit != "-" and digit > max_coord:
                temp_label[digit_index] = max_coord

        for filler in range(max_ions-num_ions):
            temp_label.extend("-"*len(sf_col))

        # Convert the label list to a string of length max_ions*len(sf_cols)
        temp_string = "".join([str(coord) for coord in temp_label])
        assert len(temp_string) == max_ions*len(sf_col)

        temp_bool = []
        for regex_id, regex in enumerate(regex_strings):
            if match(regex, temp_string) is not None:
                temp_bool.append(True)
            else:
                temp_bool.append(False)

        # Confirm there aren't multiple regex matches
        assert sum(temp_bool) < 2, \
                 "Multiple regex matches for label " +str(temp_string)

        # Here's the catch-all clause for when no regex matches
        if sum(temp_bool) == 0.0:
            temp_bool.append(True)
        else:
            temp_bool.append(False)

        # Write to filestreams if prefix is set.
        if prefix != None:
            if sum(temp_bool) == 1:
                count_files[temp_bool.index(True)].write(
                " ".join([str(col) for col in line]))
                count_files[temp_bool.index(True)].write("\n")
            else:
                count_files[len(regex_strings)].write(
                " ".join([str(col) for col in line]))
                count_files[len(regex_strings)].write("\n")

        data_output.append((temp_string, temp_bool.index(True)))

    # Close filestreams.
    if prefix != None:
        for key in count_files.keys():
            count_files[key].close()

    return data_output

# This is a script that produces a state stream from the input data
# instead of recomputing it multiple times throughout the script for things
# like transition counting. It requires sorted or sorted/padded data
# that are returned from sort_columns and pad_columns functions.
def label_columns(data_floats, num_cols=13, pad_col=4, sf_col=[5,6],
                  sort_col=3, sort_cut=0.0, max_ions=3, verbose=False):

    data_output = []
    for line in data_floats:
        temp_label = []
        num_ions = 0
        for ion in chunker(line,num_cols):

            # This functionality was implemented upon discovering that
            # zero coordinated ions can exist in either the central cavity
            # or the extracellular region of a channel. For this case,
            # a sort_value cutoff must be used to determine where the ion
            # actually is! Note, look at your SF bound ion distributions
            # to see if sort_cut is set correctly.
            sort_val = ion[sort_col]
            all_zeros = all([ion[col]==0.0 for col in sf_col])

            # I need to code a bit of logic for pre-processed data.
            if ion[pad_col] != "-":
                if (sort_val > sort_cut) and all_zeros:
                    temp_label.extend("+"*len(sf_col))
                else:
                    temp_label.extend([int(ion[col]) for col in sf_col])
            else:
                temp_label.extend("-"*len(sf_col))
            num_ions += 1

        for filler in range(max_ions-num_ions):
            temp_label.extend("-"*len(sf_col))

        # Convert the label list to a string of length max_ions*len(sf_cols)
        temp_string = "".join([str(coord) for coord in temp_label])
        data_output.append(temp_string)

    return data_output

# This script creates a uniform number of ions in the channel by paddding
# incidents where few ions are detected in the data_floats array.
def pad_columns(data_floats, num_cols=13, max_ions=3,
                time_col=0, traj_col=11, verbose=False):

    data_output = []
    for line in data_floats:

        # This is the fake line we're going to pad with.
        # We're going to assume that line[0] contains time
        # and that we have 7 additional columns other than
        # an arbitrary number of coordination count columns.
        fake_ion = [int(float(line[time_col])),0.,0.,0.] + \
                   ["-" for x in range(num_cols-7)] + \
                   [0,int(float(line[traj_col])),0]

        chunked_line = chunker(line,num_cols)
        num_ions = len(line)/num_cols

        temp_line = []
        filler_line = []

        # Note: if you have a timestep where num_ions exceeds your
        # max_ions variable, you are only going to select the
        # first "max_ions" ions and you may omit some!
        # For this reason it is best to sort your ions first,
        # so the deepest ion is preserved.
        # Make sure to run ion_counter() to see how much data loss
        # you are incurring by choosing a smaller max_ions value.
        for ion in list(chunked_line)[:max_ions]:
            temp_line.extend(ion)
            if verbose:
                print " ".join([str(ion_col) for ion_col in ion]),

        for filler in range(max_ions-num_ions):
            filler_line.extend(fake_ion)
            if verbose:
                print " ".join([str(ion_col) for ion_col in fake_ion]),

        data_output.append(temp_line+filler_line)
        if verbose:
            print "\n",

    return data_output

# This script sorts based on a particular column value, useful for
# ranking ions by z-value. Toggle sorting from highest to lowest
# with the argument plus2minus but it's pretty experimental
# when it comes to the future scripts. Namely the use of traj_col
# and going to a padded column instead of a
def sort_columns(data_floats, num_cols=13, sort_col=3,
                 verbose=False, plus2minus=True):

    data_output = []
    for line in data_floats:
        chunked_line = chunker(line,num_cols)

        temp_line = []
        temp_filler = []
        # Here we sort each ion grouping by the sort_col argument
        # with attention to reverse the list depending on plus2minus.
        for ion in sorted(chunked_line, key=lambda col: col[sort_col],
                          reverse=plus2minus):
            if ion[4] == "-":
                temp_filler.extend(ion)
            else:
                temp_line.extend(ion)

        if verbose:
            if plus2minus:
                print " ".join([str(ion_col) for ion_col in temp_line]),
                print " ".join([str(ion_col) for ion_col in temp_filler]),
            else:
                print " ".join([str(ion_col) for ion_col in temp_filler]),
                print " ".join([str(ion_col) for ion_col in temp_line]),
            print "\n",

        if plus2minus:
            data_output.append(temp_line+temp_filler)
        else:
            data_output.append(temp_filler+temp_line)

    return data_output

# This simply writes out the data_lines passed in with a simple
# ASCII format. If an filename for output is not specified it
# outputs to standard output, otherwise it outputs to the file specified.
# The write_mode argument may be used to output a series of trajectories
# to the same output file.
def write_columns(data_lines, outfile=None, write_mode="w"):
    if outfile is not None:
        fout = open(outfile,write_mode)

    for line in data_lines:
        if outfile is not None:
            fout.write(" ".join([str(col) for col in line])+"\n")
        else:
            print " ".join([str(col) for col in line])

    if outfile is not None:
        fout.close()
    return True

# This function preprocesses raw input and returns both the sorted
# ion data aswell as sorted/padded list of lists. float_cols is a list
# of column numbers that will be converted to floats, where the rest of the
# data will be converted to integers.
def process_input(filenames, sort_col=3, num_cols=13,
                  remove_frames=0, max_ions=3, traj_col=11,
                  add_time=False, time_increment=1,
                  padded=False, float_cols=[1,2,3], time_col=0):

    # There are empty lists which will contain one or multiple files
    # worth of time series data. All columns must have numeric data
    # at this stage (later you'll pad the data and add "-" strings)
    data_floats = []

    # Parse the input file and split and float all the columns.
    # This is the required format for all the functions in this
    # file.
    for filename in filenames:
        with open(filename,"r") as data_file:
            data_raw = data_file.readlines()[remove_frames:]
            data_raw_split = [line.strip().split() for line in data_raw]

        # Newer versions of my MDAnalysis code don't require adding a time
        # column (in past version only time[0] contained the timestep)
        if add_time:
            data_raw_split = add_time_column(data_raw_split, num_cols=num_cols)

        # This solves the problem of the first line not having a traj_num
        # by searching the input file, though it's a problem using a merged
        # trajectory.
        prev_traj = None
        for line in data_raw_split:
            if (len(line) > traj_col) and (len(line) > time_col):
                prev_traj = int(float(line[traj_col]))
        assert prev_traj != None, \
                 "Input file " + filename + " had no traj_id column"

        # Loop over the data and convert everything to integer
        # except the float_cols columns.
        for line_num, line in enumerate(data_raw_split):
            temp_line = []
            # This frame number will be used when the time column
            # is not located.
            frame_num = remove_frames+line_num+1
            for colindex, colvalue in enumerate(line):
                if (float_cols.count(colindex) > 0 or
                    float_cols.count(colindex % num_cols) > 0):
                    temp_line.append(float(colvalue))
                else:
                    temp_line.append(int(float(colvalue)))
            # This "pads" zero ion columns with a fake ion with
            # a timestamp.
            if len(temp_line) > 1:
                data_floats.append(temp_line)
            else:
                # TODO: There's a bug here when prev_time is
                # detected above as being half-way through a
                # file.
                data_floats.append([frame_num] +
                                   [0.,0.,0.] +
                                   ["-" for x in range(num_cols-7)] +
                                   [0,prev_traj,0])

    # TODO: Write something to remove duplicate lines.

    # Since padding normalizes your number of row entries, it's best to
    # sort your data a priori in order to capture the inner most ions
    # preferentially.
    data_floats_sorted = sort_columns(data_floats,
                                      num_cols=num_cols,
                                      sort_col=sort_col)
    if padded:
        data_floats_padded = pad_columns(data_floats_sorted,
                                     num_cols=num_cols,
                                     max_ions=max_ions,
                                     traj_col=traj_col)
        return data_floats_padded
    else:
        return data_floats_sorted

    return False

if __name__ == '__main__':
    parser = ArgumentParser(
    description='This script parses input columnular ASCII data\
    and makes it nice and pretty for subsequent analysis.')

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
    '-t', dest='traj_col', type=int, default=11,
    help='a zero inclusive column number that contains the run number')
    parser.add_argument(
    '-o', dest='outfile', type=str, default=None,
    help='the file to output the sorted padding output of all input files')
    parser.add_argument(
    '--addtime', dest='add_time', action="store_true", default=False,
    help='an optional argument to add time columns to each ion grouping')

    # The following arguments are used for regex state stream processing
    parser.add_argument(
    '-i', dest='regex', type=str, nargs="+",
    help='a list of regex values in quotes for state stream processing')
    parser.add_argument(
    '-sc', dest='sort_cut', type=float, default=0.0,
    help='a value on the sort_col range to classify zero coordinated data')
    parser.add_argument(
    '-sf', dest='sf_col', type=int, nargs="+", default=[5,6],
    help='the coordination integer columns that define the selectivity filter')
    args = parser.parse_args()

    data_f_padded = process_input(filenames=args.filenames,
                                          num_cols=args.num_cols,
                                          max_ions=args.max_ions,
                                          remove_frames=args.remove_frames,
                                          traj_col=args.traj_col,
                                          sort_col=args.sort_col,
                                          add_time=args.add_time,
                                          padded=True)

    data_f_label = label_columns(data_f_padded,
                                num_cols=args.num_cols,
                                sort_col=args.sort_col,
                                sort_cut=args.sort_cut,
                                sf_col=args.sf_col,
                                max_ions=args.max_ions)

    data_f_regex = regex_columns(data_f_padded, regex_strings=args.regex,
                                num_cols=args.num_cols,
                                sort_col=args.sort_col,
                                sort_cut=args.sort_cut,
                                sf_col=args.sf_col,
                                max_ions=args.max_ions)

    '''
    for x,y in zip(data_f_padded, data_f_regex):
        print x,y, len(y[0])
        if len(y[0]) > args.max_ions*len(args.sf_col):
            raise ValueError("State label is too long, possible bug")
    '''

    #write_columns(data_f_padded, outfile=args.outfile)

