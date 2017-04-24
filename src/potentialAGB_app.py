import numpy as np

import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt

import prepare_EOlab_layers as EO

NetCDF_file = '/home/dmilodow/DataStore_GCEL/AGB/AGBregpot.nc'
SAVEDIR = './'
dataset,geoTrans = EO.load_NetCDF(NetCDF_file)

vars = ['AGBobs','AGBpot','AGBreg']
cmaps = ['YlGnBu','YlGnBu','seismic_r']

for vv in range(0,len(vars)):
    print vars[vv]
    file_prefix = SAVEDIR + 'tropics_' + vars[vv]

    #EO.write_array_to_data_layer_GeoTiff(dataset.variables[vars[vv]],geoTrans, file_prefix)

    EO.write_array_to_display_layer_GeoTiff(dataset.variables[vars[vv]],geoTrans, file_prefix, cmap, ulim, llim,)
