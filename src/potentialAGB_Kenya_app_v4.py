"""
potentialAGB_Kenya_app_v4.py
================================================================================
Produce layers for restoration opportunity cross-comparison against other data
layers (e.g. WRI world of opportunity maps)
"""
# Import general libraries
import os
import sys
import numpy as np
import xarray as xr
import pandas as pd
import rasterio
import rasterio.mask
import fiona
# import plotting libraries
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import seaborn as sns
# import custom libraries
import prepare_EOlab_layers as EO
sys.path.append('/home/dmilodow/DataStore_DTM/FOREST2020/PotentialBiomassRFR/src')
import useful as useful
# set default cmap
plt.set_cmap(sns.light_palette('seagreen',as_cmap=True))

"""
#===============================================================================
PART A: DEFINE PATHS AND LOAD IN DATA
- Potential biomass maps (from netcdf file)
- Biome boundaries (Mapbiomas)
- WRI opportunity map
#-------------------------------------------------------------------------------
"""
country_code = 'EAFR'
country = 'Kenya'
version = '001'

path2data = '/disk/scratch/local.2/dmilodow/PotentialBiomass/processed/%s/' % country_code
path2output = '/home/dmilodow/DataStore_DTM/EOlaboratory/EOlab/KenyaPotentialAGB/'
boundaries_shp = '/home/dmilodow/DataStore_DTM/EOlaboratory/Areas/ne_50m_admin_0_tropical_countries_small_islands_removed.shp'

# load potential biomass models from netdf file
AGBpot_ds = xr.open_dataset('%s%s_%s_AGB_potential_RFR_worldclim_soilgrids_final.nc' %
                                (path2output, country_code,version))
AGBpot     = AGBpot_ds['AGBpot'].values
AGBobs     = AGBpot_ds['AGBobs'].values
AGBseq     = AGBpot-AGBobs
AGBpot_min = AGBpot_ds['AGBpot_min'].values
AGBobs_min = AGBpot_ds['AGBobs_min'].values
AGBseq_min = AGBpot_min-AGBobs_min
AGBpot_max = AGBpot_ds['AGBpot_max'].values
AGBobs_max = AGBpot_ds['AGBobs_max'].values
AGBseq_max = AGBpot_max-AGBobs_max

# Calculate cell areas
cell_areas =  useful.get_areas(latorig=AGBpot_ds.coords['lat'].values,
                            lonorig=AGBpot_ds.coords['lon'].values)
cell_areas/=10**4 # m^2 -> ha

# create national boundary mask
# - load template raster
template = rasterio.open('%s/agb/Avitabile_AGB_%s_1km.tif' % (path2data,country_code))
# - load shapefile
boundaries = fiona.open(boundaries_shp)
# - for country of interest, make mask
for feat in boundaries:
    name = feat['properties']['admin']
    if name==country:
        image,transform = rasterio.mask.mask(template,[feat['geometry']],crop=False)
masks = {}
masks[country] = image[0]>=0
masks[country] *=np.isfinite(AGBpot)

































NetCDF_file = 'kenya_ODA_PFB_mean_threshold_WorldClim2.nc'

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
dataset['seqpot_mean'][dataset['AGBpot_mean']==-9999] = -9999.

dataset['AGBpot_mean'][dataset['AGB_mean']==-9999] = -9999.
dataset['AGB_mean'][dataset['AGBpot_mean']==-9999] = -9999.

dataset['forests'][dataset['forests']!=1] = -9999.

vars = ['AGB_mean','AGBpot_mean','seqpot_mean','forests']
cmaps = ['viridis','viridis','plasma','viridis']
ulims = [50.,50.,50.,1.]
llims = [0.,0.,0.,0.]
axis_labels = ['AGB$_{obs}$ / Mg(C) ha$^{-1}$', 'AGB$_{potential}$ / Mg(C) ha$^{-1}$', 'Sequestration potential / Mg(C) ha$^{-1}$','Forest mask (1 = Forest)']


plt.figure(1, facecolor='White',figsize=[2, 1])
plt.show()


for vv in range(0,len(vars)):
    print vars[vv]
    file_prefix = SAVEDIR + 'kenya_' + vars[vv]

    # delete existing dataset if present
    if 'kenya_'+vars[vv]+'_data.tif' in os.listdir(SAVEDIR):
        os.system("rm %s" % ('kenya_'+vars[vv]+'_data.tif'))
    if 'kenya_'+vars[vv]+'_display.tif' in os.listdir(SAVEDIR):
        os.system("rm %s" % ('kenya_'+vars[vv]+'_display.tif'))

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
    file_prefix = SAVEDIR + 'kenya_' + vars[vv] + '_total'

    out_array = dataset[vars[vv]] * areas
    out_array[dataset[vars[vv]]==-9999]=-9999  # not sure why np.asarray step is needed but gets the job done
    EO.write_array_to_data_layer_GeoTiff(out_array, geoTrans, file_prefix)
    out_array=None

# Also want to write cell areas to file.  However, as this will be compared against other layers, need to carry across
# nodata values
areas_out = areas.copy()
areas_out[np.asarray(dataset[vars[0]])==-9999] = -9999
area_file_prefix = SAVEDIR + 'kenya_cell_areas'
EO.write_array_to_data_layer_GeoTiff(areas_out, geoTrans, area_file_prefix)
