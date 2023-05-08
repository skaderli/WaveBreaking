"""
This file is part of WaveBreaking.

WaveBreaking provides indices to detect, classify
and track Rossby Wave Breaking (RWB) in climate and weather data.
The tool was developed during my master thesis at the University of Bern.
Link to thesis: https://occrdata.unibe.ch/students/theses/msc/406.pdf

---
Top-level package for WaveBreaking.
"""

__author__ = "Severin Kaderli"
__license__ = "MIT"
__email__ = "severin.kaderli@unibe.ch"
__version__ = '0.3.1'

# import spatial pre-processing functions
from wavebreaking.processing.spatial import (
    calculate_momentum_flux,
    calculate_smoothed_field,
)

# import events post-processing functions
from wavebreaking.processing.events import to_xarray, track_events

# import plotting functions
from wavebreaking.processing.plots import plot_clim, plot_step, plot_tracks

# import wavebreaking indices
from wavebreaking.indices.contour_index import calculate_contours
from wavebreaking.indices.streamer_index import calculate_streamers
from wavebreaking.indices.overturning_index import calculate_overturnings
from wavebreaking.indices.cutoff_index import calculate_cutoffs
