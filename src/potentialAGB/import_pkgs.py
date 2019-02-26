# General libraries
import numpy as np
import pandas as pd
import os as os
import sys as sys

# Plotting libraries
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap as ListedColormap
import seaborn as sns
plt.cla()

# GIS libraries
import xarray as xr
import geopandas as geopandas
from geopandas.tools import sjoin
import cartopy as cartopy
from cartopy import config as config
import cartopy.crs as ccrs
from cartopy.io.shapereader import Reader as Reader
from cartopy.feature import ShapelyFeature as ShapelyFeature
from shapely.geometry import Point as Point

# custom libraries
sys.path.append('../')
import prepare_EOlab_layers as EO
sys.path.append('/home/dmilodow/DataStore_DTM/FOREST2020/EOdata/EO_data_processing/src/')
import auxilliary_functions as aux
