ChannelAnalysis: Molecular Dynamics Analysis Pipeline for Ion Channels
=======

These scripts are designed for post-processing molecular dynamics trajectories
in order to understand ion channel properties related to ion permeation. Particular
attention is paid to quantifying ionic coordination properties and permeation through
state networks of ionic binding modes. 

Functionality
=======
- ExtractProperties
  - MDAnalysis script designed to extract time-series data from multiple trajectories.
- CoordAnalysis
  - Calculates ionic occupancies in the channel.
  - Generates 1D histograms of ionic coordination numbers.
  - Generates user-defined macrostate 1D histograms.
- RotamerAnalysis (in progress)
  - Calculates rotameric populations.
  - Generates 2D histograms of Chi1-Chi2 for specified residues.
  - Generates dunking survival probabilities.

Dependencies
=======
Older versions of these packages are most likely supported, with the exception of Python 2.7
- Python 2.7
- NumPy 1.6.2
- SciPy 0.10
- MDAnalysis 0.77
- NetworkX 1.7 (Optional)
- Gnuplot (Optional)

Publication
=======
A publication which was a direct application of these scripts is forthcoming:,
- Chakrabarti, N., Ing, C., and Pomès, R., 2013 (Under Review).
