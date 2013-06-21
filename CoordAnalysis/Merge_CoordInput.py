#!/usr/bin/python
###############################################################################
#
# This script merges two MDAnalysis coordination analysis output files. In
# practice this is useful when you compute both 1st or 2nd shell coordination
# and you want to merge coordination in both of those shells. It would be
# utilized in a workflow before running Ion Preprocessor.
#
# Example: If you had multiple files with the .out extension that had
#          matching traj_id + frame_number columns, all coord_cols
#          which default to column numbers 4 5 6 7 8 9 10 would be summed:
#
#          python Merge_CoordInput.py -f *.out -o merged.out -c 13 -t 11
#
#          Note: Many important arguments are ommitted since they match
#          the default values for the test case.
#
# By Chris Ing, 2013 for Python 2.7
#
###############################################################################
from argparse import ArgumentParser
from collections import defaultdict

#a great helper function to iterate over chunks of a list
def chunker(seq, size):
    return (seq[pos:pos + size] for pos in xrange(0, len(seq), size))

# This function outputs a single file where all coord_cols are merged
# for equivalent resids at a given traj and timestep. The outputted data
# is not sorted.
def write_merged_coordination(filenames, coord_cols, outfile=None,
                              write_mode="w", traj_col=11,
                              resid_col=12, time_col=0, num_cols=13,
                              traj_digits=4, time_digits=8, remove_frames=0):
    # This datatype stores a unique traj-timeset identifier as the first key
    # and the second key is a residue id of an ion observed at that timestep
    # the value is the block of coordination and x,y,z data that are normally
    # attributed with an ion.
    merged_coord=defaultdict(dict)

    # Total file length in line numbers
    total_lines = 0

    # Parse the input file and split and float all the columns.
    # This is the required format for all the functions in this
    # file.
    for filename in filenames:
        with open(filename,"r") as data_file:
            data_raw = data_file.readlines()[remove_frames:]
            total_lines += len(data_raw)
            data_raw_split = [line.strip().split() for line in data_raw]

            # This is a traj_id detection script, we're going to pray
            # that the file input only has 1 traj_id in it or else
            # you're getting a bug.
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
                        # have more than 9999 trajectories or 99999999 frames
                        # then your data may not be correctly sorted since
                        # the identifier is a string and must be sorted
                        # alphanumerically at the end of the function.
                        u_id = ion[traj_col].split(".")[0].zfill(traj_digits)+\
                               ion[time_col].split(".")[0].zfill(time_digits)
                    elif auto_traj is not None:
                        u_id = str(auto_traj).zfill(traj_digits)+\
                               ion[time_col].split(".")[0].zfill(time_digits)
                    else:
                        raise ValueError("You input a file with no traj_id")

                    if len(line) > resid_col:
                        resid = int(float(ion[resid_col]))

                        # If that resid exists already, add the coord_cols
                        if resid in merged_coord[u_id]:
                            old_ion = merged_coord[u_id][resid]
                            for coord_col in coord_cols:
                                old_coord = float(old_ion[coord_col])
                                new_coord = float(ion[coord_col])
                                old_ion[coord_col] = str(old_coord + new_coord)
                            merged_coord[u_id][resid] = old_ion
                        # Otherwise, just create it.
                        else:
                            merged_coord[u_id][resid] = ion
                    else:
                        # In the case that no resid exists, a dictionary entry
                        # for the u_id should still be made. Though these
                        # lines will basically contribute nothing other than
                        # channel occupancy calculations.
                        merged_coord[u_id][""] = ion

    # It's entirely possible that the output may have more lines than the input
    # since we're doing a union of the two sets.
    #print len(merged_coord), total_lines

    # If a filename exists, open it, write to it, and close it.
    # Otherwise, just print.
    if outfile is not None:
        fout = open(outfile,write_mode)

    for u_id, resid_dict in sorted(merged_coord.iteritems()):
        for ions in resid_dict.values():
            if outfile == None:
                print " ".join(ions), " "
            else:
                fout.write(" ".join(ions) + " ")
        if outfile == None:
            print "\n",
        else:
            fout.write("\n")

    if outfile is not None:
        fout.close()
    return True

if __name__ == '__main__':
    parser = ArgumentParser(
    description='This script parses input columnular ASCII data\
    and makes it nice and pretty for subsequent analysis.')

    parser.add_argument(
    '-f', dest='filenames', type=str, nargs="+", required=True,
    help='a filename of coordination data from MDAnalysis trajectory data')
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
    '-coord', dest='coord_cols', type=int, nargs="+", default=[4,5,6,7,8,9,10],
    help='all the columns with coordination counts')
    parser.add_argument(
    '-o', dest='outfile', type=str, default=None,
    help='the file to output the merged input files')
    parser.add_argument(
    '-remove', dest='remove_frames', type=int, default=0,
    help='this is a number of frames to remove from the start of the data')
    args = parser.parse_args()

    write_merged_coordination(args.filenames,
                              args.coord_cols,
                              outfile=args.outfile,
                              num_cols=args.num_cols,
                              traj_col=args.traj_col,
                              time_col=args.time_col,
                              resid_col=args.resid_col,
                              remove_frames=args.remove_frames)