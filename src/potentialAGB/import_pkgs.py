# General libraries
import numpy as np
import os
import sys

# Plotting libraries
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt

# GIS libraries
import shapefile
import xarray as xr

# custom libraries
sys.path.append('../')
import prepare_EOlab_layers as EO

sys.path.append('/home/dmilodow/DataStore_DTM/FOREST2020/PotentialBiomass/src')
import geospatial_utility_tools as geo

sys.path.append('/home/dmilodow/DataStore_DTM/FOREST2020/EOdata/EO_data_processing/src/')
import auxilliary_functions as aux
