"""
potentialAGB_Brazil_app_v4.py
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
from copy import deepcopy
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
country_code = 'BRA'
country = 'Brazil'
version = '013'

path2data = '/disk/scratch/local.2/PotentialBiomass/processed/%s/' % country_code
path2model  = '/home/dmilodow/DataStore_DTM/FOREST2020/PotentialBiomassRFR/output/'
path2output = '/home/dmilodow/DataStore_DTM/EOlaboratory/EOlab/BrazilPotentialAGB/'
boundaries_shp = '/home/dmilodow/DataStore_DTM/EOlaboratory/Areas/ne_50m_admin_0_tropical_countries_small_islands_removed.shp'

source = ['globbiomass', 'avitabile']
source = ['avitabile']

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

# load opportunity map
opportunity = xr.open_rasterio('%sWRI_restoration/WRI_restoration_opportunities_%s.tif' % (path2data, country_code))[0]

# Load MapBiomas data for 2005
mb2005 = deepcopy(opportunity)
mb2005.values=useful.load_mapbiomas('BRA',timestep=20,aggregate=1)-1

for ss in source:

    # load potential biomass models from netdf file
    dataset = xr.open_dataset('%s%s_%s_AGB_potential_RFR_%s_worldclim_soilgrids_final.nc' %
                                    (path2model, country_code,version, ss))
    # Convert to C
    dataset['AGBpot'].values*=0.48
    dataset['AGBobs'].values*=0.48
    dataset['AGBpot_min'].values*=0.48
    dataset['AGBobs_min'].values*=0.48
    dataset['AGBpot_max'].values*=0.48
    dataset['AGBobs_max'].values*=0.48

    # calculate deficit aka sequestration potential
    dataset['AGBseq'] = (dataset['AGBpot']-dataset['AGBobs'])
    dataset['AGBseq_min'] = (dataset['AGBpot_min']-dataset['AGBobs_min'])
    dataset['AGBseq_max'] = (dataset['AGBpot_max']-dataset['AGBobs_max'])

    # Create potential and sequestration layers with settlements
    # maintained at original AGB (i.e. feasible restoration)
    people_mask = (mb2005.values==5)
    dataset['AGBpot_natural']=deepcopy(dataset['AGBpot'])
    dataset['AGBpot_natural'].values[people_mask]=dataset['AGBobs'].values[people_mask]
    dataset['AGBseq_natural']=deepcopy(dataset['AGBseq'])
    dataset['AGBseq_natural'].values[people_mask]=0

    """
    PART B: Create data and display layers
    - AGBobs
    - AGBpot
    - AGBseq
    - WRI restoration opportunity
    - landcover
    """
    file_prefix = path2output + country.lower() + '_'

    vars = ['AGBobs','AGBpot','AGBseq','AGBpot_natural','AGBseq_natural']
    cmaps = ['viridis','viridis','divergent','viridis','divergent']
    axis_labels = ['AGB$_{obs}$ / Mg C ha$^{-1}$', 'AGB$_{potential}$ / Mg C ha$^{-1}$', 'Sequestration potential / Mg C ha$^{-1}$', 'AGB$_{potential}$ / Mg C ha$^{-1}$', 'Sequestration potential / Mg C ha$^{-1}$']
    ulims = [200,200,100,200,200]
    llims = [0,0,-100,0,-100]
    for vv,var in enumerate(vars):
        print(var)
        if var in dataset.keys():
            file_prefix = '%s%s_%s_%s' % (path2output, country.lower(), var, ss)

            # delete existing dataset if present
            if '%s_%s_%s_data.tif' % (country.lower(),var, ss) in os.listdir(path2output):
                os.system("rm %s" % ('%s_data.tif' % (file_prefix)))
            if '%s_%s_%s_display.tif' % (country.lower(),var, ss) in os.listdir(path2output):
                os.system("rm %s" % ('%s_display.tif' % (file_prefix)))

            # apply country mask
            if ss != 'oda':
                dataset[var].values[mask==0]  = np.nan

            # write display layers
            EO.plot_legend(cmaps[vv],ulims[vv],llims[vv],axis_labels[vv], file_prefix)
            EO.write_xarray_to_display_layer_GeoTiff(dataset[vars[vv]], file_prefix, cmaps[vv], ulims[vv], llims[vv])

# WRI opportunity map
opportunity.values=opportunity.values.astype('float')
opportunity.values[mask==0]=np.nan
id =    np.arange(0,5)
labels = np.asarray( ['existing natural cover','wide-scale','mosaic','remote','urban-agriculture'])
colours = np.asarray(['#67afde', '#00883b', '#00c656', '#004c21', "#6a3b00"])

id_temp,idx_landcover,idx_id = np.intersect1d(opportunity,id,return_indices=True)
id = id[idx_id]
labels=labels[idx_id]
colours=colours[idx_id]
wri_cmap = ListedColormap(sns.color_palette(colours).as_hex())
file_prefix = '%s%s_wri' % (path2output, country.lower())
EO.plot_legend_listed(wri_cmap,labels,'',file_prefix,figsize=[2,1])

if '%s_wri_data.tif' % (country.lower()) in os.listdir(path2output):
    os.system("rm %s" % ('%s_data.tif' % (file_prefix)))
if '%s_wri_display.tif' % (country.lower()) in os.listdir(path2output):
    os.system("rm %s" % ('%s_display.tif' % (file_prefix)))
EO.write_xarray_to_display_layer_GeoTiff(opportunity, file_prefix, wri_cmap, 4, 0)

# Mapbiomas land cover data
lc_class = np.array(['Natural Forest','Natural Non-Forest','Plantation','Pasture','Agriculture','Urban','Other'])
colours = np.asarray(['#1f4423', '#bbfcac', '#935132', '#ffd966', '#e974ed','#af2a2a','#d5d5e5'])
lc_id = np.arange(0,7)

id_temp,idx_landcover,idx_id = np.intersect1d(mb2005,lc_id,return_indices=True)
lc_id = lc_id[idx_id]
lc_class=lc_class[idx_id]
colours=colours[idx_id]
mb_cmap = ListedColormap(sns.color_palette(colours).as_hex())
mb_cmap_rev = ListedColormap(sns.color_palette(colours[::-1]).as_hex())
file_prefix = '%s%s_mapbiomas' % (path2output, country.lower())
EO.plot_legend_listed(mb_cmap_rev,lc_class[::-1],'',file_prefix,figsize=[2,2])

file_prefix = '%s%s_mapbiomas_lc_2005' % (path2output, country.lower())
if '%s_mapbiomas_lc_2005_data.tif' % (country.lower()) in os.listdir(path2output):
    os.system("rm %s" % ('%s%s_data.tif' % (file_prefix,country.lower())))

if '%s_mapbiomas_lc_2005_display.tif' % (country.lower()) in os.listdir(path2output):
    os.system("rm %s" % ('%s%s_display.tif' % (file_prefix,country.lower())))
EO.write_xarray_to_display_layer_GeoTiff(mb2005, file_prefix, mb_cmap, 6, 0)
