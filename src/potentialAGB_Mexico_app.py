import numpy as np
import os
import sys

import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt

import prepare_EOlab_layers as EO

sys.path.append('/home/dmilodow/DataStore_DTM/FOREST2020/PotentialBiomass/src')
import geospatial_utility_tools as geo


# Get perceptionally uniform colourmaps
sys.path.append('/home/dmilodow/DataStore_DTM/FOREST2020/EOdata/EO_data_processing/src/plot_EO_data/colormap/')
import colormaps as cmaps
plt.register_cmap(name='viridis', cmap=cmaps.viridis)
plt.register_cmap(name='inferno', cmap=cmaps.inferno)
plt.register_cmap(name='plasma', cmap=cmaps.plasma)
plt.register_cmap(name='magma', cmap=cmaps.magma)
plt.set_cmap(cmaps.viridis)

DATADIR = '/disk/scratch/local.2/mexico_PFB/'
SAVEDIR = '/home/dmilodow/DataStore_DTM/EOlaboratory/EOlab/MexicoPotentialAGB/'
NetCDF_file = 'mexico_PFB.nc'

ds,geoTrans = EO.load_NetCDF(DATADIR+NetCDF_file,lat_var = 'lat', lon_var = 'lon')
resampling_scalar = 3.
vars = ['AGB-MEX','AGBpot','forests']
dataset, geoTrans = EO.resample_dataset(ds,geoTrans,vars,resampling_scalar)

dataset['AGBreg'] = dataset['AGBpot']-dataset['AGB-MEX']
dataset['AGBreg'][dataset['AGB-MEX']==-9999] = -9999.

dataset['refpot'] = dataset['AGBpot']-dataset['AGB-MEX']
dataset['refpot'][dataset['AGB-MEX']==-9999] = -9999.
dataset['refpot'][dataset['forests']==1] = 0.

vars = ['AGB-MEX','AGBpot','AGBreg','refpot','forests']
cmaps = ['viridis','viridis','PRGn','inferno','viridis']
ulims = [230.,230.,100.,100.,1.]
llims = [0.,0.,-100.,0.,0.]
axis_labels = ['AGB$_{obs}$ / Mg(C) ha$^{-1}$', 'AGB$_{potential}$ / Mg(C) ha$^{-1}$', 'AGB deficit / Mg(C) ha$^{-1}$', 'Reforestation potential / Mg(C) ha$^{-1}$','Forest mask (1 = Forest)']

for vv in range(0,len(vars)):
    print vars[vv]
    file_prefix = SAVEDIR + 'mexico_' + vars[vv]

    # delete existing dataset if present
    if 'mexico_'+vars[vv]+'_data.tif' in os.listdir(SAVEDIR):
        os.system("rm %s" % ('mexico_'+vars[vv]+'_data.tif'))
    if 'mexico_'+vars[vv]+'_display.tif' in os.listdir(SAVEDIR):
        os.system("rm %s" % ('mexico_'+vars[vv]+'_display.tif'))

    EO.write_array_to_display_layer_GeoTiff(dataset[vars[vv]], geoTrans, file_prefix, cmaps[vv], ulims[vv], llims[vv])
    EO.plot_legend(cmaps[vv],ulims[vv],llims[vv],axis_labels[vv], file_prefix)

# Now write the quantitative display layers to file.  These probide the 
rows, cols = dataset[vars[0]].shape
latitude = np.arange(geoTrans[3],rows*geoTrans[5]+geoTrans[3],geoTrans[5])
longitude =  np.arange(geoTrans[0],cols*geoTrans[1]+geoTrans[0],geoTrans[1])
areas = geo.calculate_cell_area_array(latitude,longitude, area_scalar = 1./10.**4,cell_centred=False)

# loop through the variables, multiplying by cell areas to give values in Mg
for vv in range(0,len(vars)):
    print vars[vv]
    file_prefix = SAVEDIR + 'mexico_' + vars[vv] + '_total_data'

    out_array = dataset[vars[vv]] * areas
    out_array[dataset[vars[vv]]==-9999]=-9999  # not sure why np.asarray step is needed but gets the job done
    EO.write_array_to_data_layer_GeoTiff(out_array, geoTrans, file_prefix)
    out_array=None

# Also want to write cell areas to file.  However, as this will be compared against other layers, need to carry across
# nodata values
areas_out = areas.copy()
areas_out[np.asarray(dataset[vars[0]])==-9999] = -9999
area_file_prefix = SAVEDIR + 'mexico_cell_areas_data'
EO.write_array_to_data_layer_GeoTiff(areas_out, geoTrans, area_file_prefix)
