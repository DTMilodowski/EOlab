import numpy as np
import os
import sys

import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt

import prepare_EOlab_layers as EO

sys.path.append('/home/dmilodow/DataStore_DTM/FOREST2020/PotentialBiomass/src')
import geospatial_utility_tools as geo

DATADIR = '/disk/scratch/local.2/kenya_PFB/'
SAVEDIR = '/home/dmilodow/DataStore_DTM/EOlaboratory/EOlab/KenyaPotentialAGB/'
site = 'kenya'
ds,geoTrans = EO.load_NetCDF(NetCDF_file,lat_var = 'lat', lon_var = 'lon')
resampling_scalar = 3.

#vars = ['AGBobs','AGBpot','AGBreg']

dataset, geoTrans = EO.resample_dataset(ds,geoTrans,vars,resampling_scalar)

cmaps = ['viridis','viridis','PRGn']
ulims = [230.,230.,100.]
llims = [0.,0.,-100.]
axis_labels = ['AGB$_{obs}$ / Mg(C) ha$^{-1}$', 'AGB$_{potential}$ / Mg(C) ha$^{-1}$', 'AGB deficit / Mg(C) ha$^{-1}$']

for vv in range(0,len(vars)):
    print vars[vv]
    file_prefix = SAVEDIR + site + '_' + vars[vv]

    # delete existing dataset if present
    if site+'_'+vars[vv]+'_data.tif' in os.listdir(SAVEDIR):
        os.system("rm %s" % (site+'_'+vars[vv]+'_data.tif'))
    if site+'_'+vars[vv]+'_display.tif' in os.listdir(SAVEDIR):
        os.system("rm %s" % (site+'_'+vars[vv]+'_display.tif'))

    #EO.write_array_to_data_layer_GeoTiff(dataset.variables[vars[vv]],geoTrans, file_prefix)

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
    file_prefix = SAVEDIR + site + '_' + vars[vv] + '_total_data'

    out_array = dataset[vars[vv]] * areas
    out_array[dataset[vars[vv]]==-9999]=-9999  # not sure why np.asarray step is needed but gets the job done
    EO.write_array_to_data_layer_GeoTiff(out_array, geoTrans, file_prefix)
    out_array=None

# Also want to write cell areas to file.  However, as this will be compared against other layers, need to carry across
# nodata values
areas_out = areas.copy()
areas_out[np.asarray(dataset[vars[0]])==-9999] = -9999
area_file_prefix = SAVEDIR + site
EO.write_array_to_data_layer_GeoTiff(areas_out, geoTrans, area_file_prefix)
