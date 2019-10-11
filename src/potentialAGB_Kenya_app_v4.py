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
from matplotlib.colors import ListedColormap
import seaborn as sns
# import custom libraries
import prepare_EOlab_layers as EO
sys.path.append('/home/dmilodow/DataStore_DTM/FOREST2020/PotentialBiomassRFR/src')
import useful as useful
# set default cmap
plt.set_cmap('viridis')
plt.register_cmap(name='divergent', cmap=sns.diverging_palette(275,150,l=66,s=90,as_cmap=True))
#sns.light_palette('seagreen',as_cmap=True)
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
path2model  = '/home/dmilodow/DataStore_DTM/FOREST2020/PotentialBiomassRFR/output/'
path2output = '/home/dmilodow/DataStore_DTM/EOlaboratory/EOlab/KenyaPotentialAGB/'
boundaries_shp = '/home/dmilodow/DataStore_DTM/EOlaboratory/Areas/ne_50m_admin_0_tropical_countries_small_islands_removed.shp'

# load potential biomass models from netdf file
dataset = xr.open_dataset('%s%s_%s_AGB_potential_RFR_worldclim_soilgrids_final.nc' %
                                (path2model, country_code,version))
dataset['AGBseq'] = (dataset['AGBpot']-dataset['AGBobs'])/2.
dataset['AGBseq_min'] = (dataset['AGBpot_min']-dataset['AGBobs_min'])/2.
dataset['AGBseq_max'] = (dataset['AGBpot_max']-dataset['AGBobs_max'])/2.

# load opportunity map
opportunity = xr.open_rasterio('%sWRI_restoration/WRI_restoration_opportunities_%s.tif' % (path2data, country_code))[0]

# Load ESACCI data for 2005
esacci2005 = useful.load_esacci('EAFR',year=2005,aggregate=1)

# Create potential and sequestration layers with agriculture and settlements
# maintained at original AGB (i.e. feasible restoration)
people_mask = np.any((esacci2005==4,esacci2005==1),axis=0)
dataset['AGBpot_natural']=dataset['AGBpot'].copy()
dataset['AGBpot_natural'].values[people_mask]=dataset['AGBobs'].values[people_mask]
dataset['AGBseq_natural']=dataset['AGBseq'].copy()
dataset['AGBseq_natural'].values[people_mask]=0


# create and apply national boundary mask
# - load template raster
template = rasterio.open('%s/agb/Avitabile_AGB_%s_1km.tif' % (path2data,country_code))
# - load shapefile
boundaries = fiona.open(boundaries_shp)
# - for country of interest, make mask
mask = np.zeros(template.shape)
for feat in boundaries:
    name = feat['properties']['admin']
    if name==country:
        image,transform = rasterio.mask.mask(template,[feat['geometry']],crop=False)
        mask[image[0]>=0]=1

"""
PART B: Create data and display layers
- AGBobs
- AGBpot
- AGBseq
- WRI restoration opportunity
- ?landcover
"""
file_prefix = path2output + country.lower() + '_'

vars = ['AGBobs','AGBpot','AGBseq','AGBpot_natural','AGBseq_natural']
cmaps = ['viridis','viridis','divergent','viridis','divergent']
axis_labels = ['AGB$_{obs}$ / Mg ha$^{-1}$', 'AGB$_{potential}$ / Mg ha$^{-1}$', 'Sequestration potential / Mg(C) ha$^{-1}$', 'AGB$_{potential}$ / Mg ha$^{-1}$', 'Sequestration potential / Mg(C) ha$^{-1}$'']
ulims = [300,300,150,300,150]
llims = [0,0,-150,0,-150]
for vv,var in enumerate(vars):
    print(var)
    file_prefix = '%s%s_%s' % (path2output, country.lower(), var)

    # delete existing dataset if present
    if '%s_%s_data.tif' % (country.lower(),var) in os.listdir(path2output):
        os.system("rm %s" % ('%s_data.tif' % (file_prefix)))
    if '%s_%s_display.tif' % (country.lower(),var) in os.listdir(path2output):
        os.system("rm %s" % ('%s_display.tif' % (file_prefix)))

    # apply country mask
    dataset[var].values[mask==0]  = np.nan

    # write display layers
    EO.plot_legend(cmaps[vv],ulims[vv],llims[vv],axis_labels[vv], file_prefix)
    EO.write_xarray_to_display_layer_GeoTiff(dataset[vars[vv]], file_prefix, cmaps[vv], ulims[vv], llims[vv])


# WRI opportunity map
opportunity.values=opportunity.values.astype('float')
opportunity.values[mask==0]=np.nan
id =    np.arange(0,5)
#labels = np.asarray( ['existing forest','wide-scale','mosaic','remote','urban-agriculture'])
#colours = np.asarray(['#00ac4b', '#055f9a', '#318ac4', '#023a5f', "#955300"])
labels = np.asarray( ['existing natural cover','wide-scale','mosaic','remote','urban-agriculture'])
colours = np.asarray(['#67afde', '#00883b', '#00c656', '#004c21', "#6a3b00"])

id_temp,idx_landcover,idx_id = np.intersect1d(opportunity,id,return_indices=True)
id = id[idx_id]
labels=labels[idx_id]
colours=colours[idx_id]
wri_cmap = ListedColormap(sns.color_palette(colours).as_hex())
EO.plot_legend_listed(wri_cmap,labels,'',file_prefix,figsize=[2,1])

file_prefix = '%s%s_wri' % (path2output, country.lower())
if '%s_wri_data.tif' % (country.lower()) in os.listdir(path2output):
    os.system("rm %s" % ('%s%s_data.tif' % (file_prefix,country.lower())))

if '%s_wri_display.tif' % (country.lower()) in os.listdir(path2output):
    os.system("rm %s" % ('%s%s_display.tif' % (path2output,country.lower())))
EO.write_xarray_to_display_layer_GeoTiff(opportunity, file_prefix, wri_cmap, 4, 0)

# ESA CCI Land cover
esacci2005.values=esacci2005.values.astype('float')
esacci2005.values[mask==0]=np.nan
lc = essacci2005.values.copy()
esacci2005.values[lc==2] = 0
esacci2005.values[lc==3] = 1
esacci2005.values[lc==6] = 2
esacci2015.values[lc==8] = 3
esacci2015.values[lc==9] = 4
esacci2015.values[lc==4] = 5
esacci2015.values[lc==1] = 6
lc_class = ['Forest','Grass','Shrub','Sparse','Bare','Wetland','Agriculture','Urban']
colours = np.asarray(['#67afde', '#00883b', '#00c656', '#004c21', "#6a3b00"])
lc_id = np.arange(0,7)

id_temp,idx_landcover,idx_id = np.intersect1d(esacci2005.values,id,return_indices=True)
id = id[idx_id]
lc_class=lc_class[idx_id]
colours=colours[idx_id]
esacci_cmap = ListedColormap(sns.color_palette(colours).as_hex())
EO.plot_legend_listed(esacci_cmap,labels,'',file_prefix,figsize=[2,1])

file_prefix = '%s%s_esacci_lc_2005' % (path2output, country.lower())
if '%s_wri_data.tif' % (country.lower()) in os.listdir(path2output):
    os.system("rm %s" % ('%s%s_data.tif' % (file_prefix,country.lower())))

if '%s_wri_display.tif' % (country.lower()) in os.listdir(path2output):
    os.system("rm %s" % ('%s%s_display.tif' % (path2output,country.lower())))
EO.write_xarray_to_display_layer_GeoTiff(esacci2005, file_prefix, esacci2005_cmap, 6, 0)
