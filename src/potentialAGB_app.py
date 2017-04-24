import numpy as np

import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt

import prepare_EOlab_layers as EO

NetCDF_file = '/home/dmilodow/DataStore_GCEL/AGB/AGBregpot.nc'
dataset,geoTrans = EO.load_NetCDF(NetCDF_file)

vars = ['AGBobs','AGBpot','AGBreg']

