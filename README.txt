ChannelAnalysis: Molecular Dynamics Analysis Pipeline for Ion Channels
=======

These scripts are designed for post-processing molecular dynamics trajectories
in order to understand ion channel properties related to ion permeation. Particular
attention is paid to quantifying ionic coordination properties and permeation through
state networks of ionic binding modes.

Functionality
=======
- ExtractProperties (Forthcoming)
  - MDAnalysis script designed to extract time-series data from multiple trajectories for CoordAnalysis/RotamerAnalysis.
- CoordAnalysis
  - Produces timeseries plots of ionic positions over time stacked with coordination data.
  - Generates 1D histograms of channel ion occupancy states along with populations.
  - Generates 1D histograms of ionic coordination numbers.
  - Generates 1D histograms of ionic binding modes as a function of axial position.
  - Generates 1D histograms of regular-expression defined macrostates.
  - Generates 2D histograms for axial positions of ion pairs, as well as 1D histograms of pair distance.
  - Calculates microstate to microstate transitions as well as intermediate-mediated transitions
  - Generates Gephi-readable graph files for all microstates.
  - Generates network graph for macrostate transitions.
  - Extra Utility to merge coordinate input (1st and 2nd shell coordination).
- RotamerAnalysis
  - Calculates rotameric populations.
  - Generates 1D histograms of Chi2 for specified residues.
  - Generates 2D histograms of Chi1-Chi2 for specified residues.
  - Generates rotamer state survival probabilities and fits to exponential decay functions.
- PoreAnalysis
  - Produces histogram of channel atom positions.

Future Functionality
=======
- PoreAnalysis scripts that track pore dynamics for correlation to ionic coordination.
- Configuration file for all input flags (so they are saved outside the command-line).
- Matplotlib plotting for all figures to create "infocards" on each dataset.
- iPython Notebook Support for "Init" scripts.
- Test dataset for working example on Wakari.

Dependencies
=======
Older versions of these packages are most likely supported, with the exception of Python 2.7
- Python 2.7
- NumPy 1.6.2
- SciPy 0.10
- MDAnalysis 0.77
- Matplotlib 1.2.1 (Optional for plotting 1D histograms)
- NetworkX 1.7 (Optional for microstate graphs)
- Gnuplot (Optional for plotting 1D histograms)

Publication
=======
A publication with direct application of these scripts is forthcoming:,
- Chakrabarti, N., Ing, C., Payandeh, J., Zheng, N., Catterall, W.A., and Pomès, R., PNAS, 2013 (In Press).
