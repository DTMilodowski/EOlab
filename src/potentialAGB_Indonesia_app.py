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

plt.figure(1, facecolor='White',figsize=[2, 1])
plt.show()

DATADIR = '/disk/scratch/local.2/southeast_asia_PFB/'
SAVEDIR = '/home/dmilodow/DataStore_DTM/EOlaboratory/EOlab/IndonesiaPotentialAGB/v1.0/'
NetCDF_file = 'southeast_asia_PFB_mean.nc'

ds,geoTrans = EO.load_NetCDF(DATADIR+NetCDF_file,lat_var = 'lat', lon_var = 'lon')
resampling_scalar = 3.
vars = ['AGB_mean','AGBpot_mean','forests']
dataset, geoTrans = EO.resample_dataset(ds,geoTrans,vars,resampling_scalar)

# sequestration potential is defined by pixels with positive potential biomass that
# are not already forests
dataset['seqpot_mean'] = dataset['AGBpot_mean']-dataset['AGB_mean']
dataset['seqpot_mean'][dataset['forests']==1] = 0.
dataset['seqpot_mean'][dataset['seqpot_mean']<0] = 0.
dataset['seqpot_mean'][dataset['AGB_mean']==-9999] = -9999.

dataset['forests'][dataset['forests']!=1] = -9999.

vars = ['AGB_mean','AGBpot_mean','seqpot_mean','forests']
cmaps = ['viridis','viridis','plasma','viridis']
ulims = [200.,200.,100.,1.]
llims = [0.,0.,0.,0.]
axis_labels = ['AGB$_{obs}$ / Mg(C) ha$^{-1}$', 'AGB$_{potential}$ / Mg(C) ha$^{-1}$', 'Sequestration potential / Mg(C) ha$^{-1}$', 'Forest mask (1 = Forest)']

for vv in range(0,len(vars)):
    print vars[vv]
    file_prefix = SAVEDIR + 'se_asia_' + vars[vv]

    # delete existing dataset if present
    if 'se_asia_'+vars[vv]+'_data.tif' in os.listdir(SAVEDIR):
        os.system("rm %s" % (SAVEDIR+'se_asia_'+vars[vv]+'_data.tif'))
    if 'se_asia_'+vars[vv]+'_display.tif' in os.listdir(SAVEDIR):
        os.system("rm %s" % (SAVEDIR+'se_asia_'+vars[vv]+'_display.tif'))

    EO.write_array_to_display_layer_GeoTiff(dataset[vars[vv]], geoTrans, file_prefix, cmaps[vv], ulims[vv], llims[vv])
    EO.plot_legend(cmaps[vv],ulims[vv],llims[vv],axis_labels[vv], file_prefix)

rows, cols = dataset[vars[0]].shape
latitude = np.arange(geoTrans[3],rows*geoTrans[5]+geoTrans[3],geoTrans[5])
longitude =  np.arange(geoTrans[0],cols*geoTrans[1]+geoTrans[0],geoTrans[1])
areas = geo.calculate_cell_area_array(latitude,longitude, area_scalar = 1./10.**4,cell_centred=False)

# loop through the variables, multiplying by cell areas to give values in Mg
for vv in range(0,len(vars)):
    print vars[vv]
    
    if 'se_asia_'+vars[vv]+'_total_data.tif' in os.listdir(SAVEDIR):
        os.system("rm %s" % (SAVEDIR+'se_asia_'+vars[vv]+'_total_data.tif'))
        
    file_prefix = SAVEDIR + 'se_asia_' + vars[vv] + '_total'

    out_array = dataset[vars[vv]] * areas
    out_array[dataset[vars[vv]]==-9999]=-9999
    EO.write_array_to_data_layer_GeoTiff(out_array, geoTrans, file_prefix)
    out_array=None

# Also want to write cell areas to file.  However, as this will be compared against other layers, need to carry across
# nodata values
areas_out = areas.copy()
areas_out[np.asarray(dataset[vars[0]])==-9999] = -9999
if 'se_asia_cell_areas_data.tif' in os.listdir(SAVEDIR):
    os.system("rm %s" % (SAVEDIR+'se_asia_cell_areas_data.tif'))
area_file_prefix = SAVEDIR + 'se_asia_cell_areas'
EO.write_array_to_data_layer_GeoTiff(areas_out, geoTrans, area_file_prefix)



"""
# Finally, we also want to write quantitative display layers for the maximum and minimum biomass,
# potential biomass and sequestration potential so that we can define uncertainy boundaries.
NetCDF_file = 'indonesia_PFB_lower.nc'
ds,geoTrans = EO.load_NetCDF(DATADIR+NetCDF_file,lat_var = 'lat', lon_var = 'lon')
resampling_scalar = 3.
vars = ['AGB_lower','AGBpot_lower','forests']
dataset, geoTrans = EO.resample_dataset(ds,geoTrans,vars,resampling_scalar)

# sequestration potential is defined by pixels with positive potential biomass that
# are not already forests
dataset['seqpot_lower'] = dataset['AGBpot_lower']-dataset['AGB_lower']
dataset['seqpot_lower'][dataset['AGB_lower']==-9999] = -9999.
dataset['seqpot_lower'][dataset['forests']==1] = 0.
dataset['seqpot_lower'][dataset['seqpot_lower']<0] = 0.

dataset['forests'][dataset['forests']!=1] = -9999.

vars = ['AGB_lower','AGBpot_lower','seqpot_lower']
for vv in range(0,len(vars)):
    print vars[vv]
    
    if 'indonesia_'+vars[vv]+'_total_data.tif' in os.listdir(SAVEDIR):
        os.system("rm %s" % (SAVEDIR+'indonesia_'+vars[vv]+'_total_data.tif'))
        
    file_prefix = SAVEDIR + 'indonesia_' + vars[vv] + '_total_data'

    out_array = dataset[vars[vv]] * areas
    out_array[dataset[vars[vv]]==-9999]=-9999 
    EO.write_array_to_data_layer_GeoTiff(out_array, geoTrans, file_prefix)
    out_array=None

NetCDF_file = 'indonesia_PFB_upper.nc'
ds,geoTrans = EO.load_NetCDF(DATADIR+NetCDF_file,lat_var = 'lat', lon_var = 'lon')
resampling_scalar = 3.
vars = ['AGB_upper','AGBpot_upper','forests']
dataset, geoTrans = EO.resample_dataset(ds,geoTrans,vars,resampling_scalar)

# sequestration potential is defined by pixels with positive potential biomass that
# are not already forests
dataset['seqpot_upper'] = dataset['AGBpot_upper']-dataset['AGB_upper']
dataset['seqpot_upper'][dataset['AGB_upper']==-9999] = -9999.
dataset['seqpot_upper'][dataset['forests']==1] = 0.
dataset['seqpot_upper'][dataset['seqpot_upper']<0] = 0.

dataset['forests'][dataset['forests']!=1] = -9999.

vars = ['AGB_upper','AGBpot_upper','seqpot_upper']
for vv in range(0,len(vars)):
    print vars[vv]
    
    if 'indonesia_'+vars[vv]+'_total_data.tif' in os.listdir(SAVEDIR):
        os.system("rm %s" % (SAVEDIR+'indonesia_'+vars[vv]+'_total_data.tif'))
        
    file_prefix = SAVEDIR + 'indonesia_' + vars[vv] + '_total_data'

    out_array = dataset[vars[vv]] * areas
    out_array[dataset[vars[vv]]==-9999]=-9999  
    EO.write_array_to_data_layer_GeoTiff(out_array, geoTrans, file_prefix)
    out_array=None
"""
