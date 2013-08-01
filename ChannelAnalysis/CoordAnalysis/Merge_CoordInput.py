#!/usr/bin/python
###############################################################################
#
# This script merges two MDAnalysis coordination analysis output files. In
# practice this is useful when you have data from two ion species. It would be
# utilized in a workflow before running Ion Preprocessor to seperate
# data into distinct groups.
#
# Example: If you had multiple files with the .out extension that had
#          matching traj_id + frame_number columns, atoms would be appended
#          at that timestep.
#
#          python Merge_CoordInput.py -f *.out -o merged.out -c 13 -t 11
#
# By Chris Ing, 2013 for Python 2.7
#
###############################################################################
from argparse import ArgumentParser
from collections import defaultdict

#a great helper function to iterate over chunks of a list
def chunker(seq, size):
    return (seq[pos:pos + size] for pos in xrange(0, len(seq), size))

# This function outputs trajectory files for each set of traj_ids
# with the ions appended to the same line. The resid column now encodes
# which block of datafiles the ion came from with a trailing integer (0 or 1)
def write_merged_coordination(filenames0, filenames1, prefix=None,
                              write_mode="w", traj_col=11,
                              resid_col=12, time_col=0, num_cols=13,
                              traj_digits=4, time_digits=8, remove_frames=0):
    # This datatype stores a unique traj-timeset identifier as the first key
    # and the second key is a residue id of an ion observed at that timestep
    # the value is the block of coordination and x,y,z data that are normally
    # attributed with an ion.
    merged_coord=defaultdict(list)

    # Total file length in line numbers
    total_lines = 0

    # As you can see, I've hardcoded two sets of filenames that correspond
    # to ionic coordination for two ion species.
    for species_num, filenames in enumerate([filenames0, filenames1]):
        for filename in filenames:
            with open(filename,"r") as data_file:
                data_raw = data_file.readlines()[remove_frames:]
                total_lines += len(data_raw)
                data_raw_split = [line.strip().split() for line in data_raw]

                # This is a traj_id detection script, we're going to pray
                # that the file input only has at least 1 traj_id in it or else
                # you're getting a bug. Note that if a whole trajectory has no
                # ion occupancy than this is the only issue.
                auto_traj = None
                for line in data_raw_split:
                    if (len(line) > traj_col) and (len(line) > time_col):
                        auto_traj = int(float(line[traj_col]))

                # Now loop over the data and convert everything to integer
                # except the float_cols columns.
                for line in data_raw_split:
                    # Loop over all ion's at that timestep
                    for ion in chunker(line,num_cols):
                        # The default behavior is to read the column traj_id
                        # but in the case that it doesn't exist, use the auto
                        #
                        if (len(line) > traj_col) and (len(line) > time_col):
                            # Prepare a unique identifier for that line, if you
                            # have more than 9999 trajectories or 99999999
                            # frames then your data may not be correctly sorted
                            # since the identifier is a string and must be
                            # sorted alphanumerically at the end of the
                            # function.
                            u_id_traj = ion[traj_col].split(".")[0]
                            u_id_time = ion[time_col].split(".")[0]
                            u_id = u_id_traj.zfill(traj_digits) + "-" + \
                                   u_id_time.zfill(time_digits)
                        elif auto_traj is not None:
                            u_id_traj = str(auto_traj)
                            u_id_time = ion[time_col].split(".")[0]
                            u_id = u_id_traj.zfill(traj_digits) + "-" + \
                                   u_id_time.zfill(time_digits)
                        else:
                            raise ValueError("input file detected no traj_id")

                        ion[resid_col] = ion[resid_col]+str(species_num)
                        merged_coord[u_id].append(ion)

    file_dict = {}
    for u_id, resid_list in sorted(merged_coord.iteritems()):

        if prefix is not None:
            traj_id = u_id.split("-")[0]
            if traj_id not in file_dict:
                file_dict[traj_id] = open(prefix+"_n"+traj_id,write_mode)

        for ions in resid_list:
            if prefix is not None:
                file_dict[traj_id].write(" ".join(ions) + " ")
            else:
                print " ".join(ions), " "

        if prefix is not None:
            file_dict[traj_id].write("\n")
        else:
            print "\n",

    if prefix is not None:
        file_dict[traj_id].close()

    return True

if __name__ == '__main__':
    parser = ArgumentParser(
    description='This script parses input columnular ASCII data\
    and makes it nice and pretty for subsequent analysis.')

    parser.add_argument(
    '-f1', dest='filenames1', type=str, nargs="+", required=True,
    help='filenames of coordination data for ion species 1')
    parser.add_argument(
    '-f2', dest='filenames2', type=str, nargs="+", required=True,
    help='filenames of coordination data for ion species 2')
    parser.add_argument(
    '-c', dest='num_cols', type=int, default=13,
    help='the number of columns per ion in the input')
    parser.add_argument(
    '-r', dest='resid_col', type=int, default=12,
    help='a zero inclusive column number that contains the ion resid')
    parser.add_argument(
    '-t', dest='traj_col', type=int, default=11,
    help='a zero inclusive column number that contains the run number')
    parser.add_argument(
    '-time', dest='time_col', type=int, default=0,
    help='a zero inclusive column number that contains the frame number')
    parser.add_argument(
    '-p', dest='prefix', type=str, default="merged_ion_species",
    help='this is the prefix of the histogram output')
    parser.add_argument(
    '-remove', dest='remove_frames', type=int, default=0,
    help='this is a number of frames to remove from the start of the data')
    args = parser.parse_args()

    write_merged_coordination(args.filenames1, args.filenames2,
                              prefix=args.prefix,
                              num_cols=args.num_cols,
                              traj_col=args.traj_col,
                              time_col=args.time_col,
                              resid_col=args.resid_col,
                              remove_frames=args.remove_frames)