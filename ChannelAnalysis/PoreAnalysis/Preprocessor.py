#!/usr/bin/python
###############################################################################
#
# Accepts multiple files of atomic positions of pore lining residues and
# processes them for plotting and producing correlations to other
# ChannelAnalysis time series data.
#
# Example: For 14-column data with this type (described elsewhere):
#
# 1 7.0 0.413 0.373 0.294 0.300 0.282 0.425 0.358 0.246 0.422 0.305 0.392 0.350
# 2 7.0 0.412 0.337 0.280 0.388 0.292 0.419 0.384 0.233 0.469 0.287 0.389 0.301
#
#          The following command will load that datafile into memory, strip
#          the first 2000 lines and ready that data for subsequent analysis.
#
#          python ChannelAtom_Preprocessor.py -f nav.n7.thr -remove 2000
#
# This script is typically used during figure plotting to obtain histograms
# of channel lining oxygens.
#
# By Chris Ing, 2013 for Python 2.7
#
###############################################################################
from argparse import ArgumentParser
import gzip

# A helper function to open with gzip if the file is gzipped
def file_opener(fname):
    if fname.endswith('.gz'):
        return gzip.open(fname)
    else:
        return open(fname)

def process_channelatoms(filenames, remove_frames=0):

    # There are empty lists which will contain one or multiple files
    # worth of time series data. All columns must have numeric data
    # at this stage (later you'll pad the data and add "-" strings)
    data_floats = []

   # Parse the input file and split and float all the columns.
    # This is the required format for all the functions in this
    # file.
    for filename in filenames:
        with file_opener(filename) as data_file:
            data_raw = data_file.readlines()[remove_frames:]
            data_raw_split = [line.strip().split() for line in data_raw]

            for line in data_raw_split:
                data_floats.append([float(col) for col in line])

    return data_floats

if __name__ == '__main__':
    parser = ArgumentParser(
    description='This script parses input columnular ASCII data\
    of channel atoms and makes it nice and pretty for subsequent analysis.')

    parser.add_argument(
    '-f', dest='filenames', type=str, nargs="+", required=True,
    help='a filename of atom data from MDAnalysis trajectory data')
    parser.add_argument(
    '-remove', dest='remove_frames', type=int, default=0,
    help='this is a number of frames to remove from the start of the data')
    args = parser.parse_args()

    data_f_processed = process_channelatoms(filenames=args.filenames,
                                            remove_frames=args.remove_frames)

    print data_f_processed