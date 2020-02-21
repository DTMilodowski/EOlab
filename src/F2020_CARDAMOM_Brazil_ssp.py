"""
F2020_CARDAMOM_Brazil_ssp.py
================================================================================
This contains code to produce EOlab layers for Brazil Forests2020 applications
that visualise key outcomes from CARDAMOM simulations that forecast the carbon
cycle dynamics under different reference scenarios.
Output display layers include:
1) Mean initial values (wood, soil)
2) Mean 2050 values (wood, soil)
3) Change in AGB by 2050
4) Change in AGB (greyscale if significance level not sufficient)
Significance will in this case be specified based on a pixel mask that has been
calculated prior to running this script.
--------------------------------------------------------------------------------
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
import colour_tools as clt
sys.path.append('/home/dmilodow/DataStore_DTM/FOREST2020/PotentialBiomassRFR/src')
import useful as useful
# set default cmap

cm1 = sns.cubehelix_palette(256,start=2.33,rot=-0.33,dark=0.2,light=0.95,as_cmap=True)
cm1.name='foliage'
plt.register_cmap(name=cm1.name, cmap=cm1)
rotation = 0.5
cm2=sns.cubehelix_palette(128,start=0.5,rot=0.5/0.95,dark=0.2,light=0.95,reverse=True,as_cmap=True)
cm3=sns.cubehelix_palette(128,start=0.5,rot=-0.5/0.95,dark=0.2,light=0.95,reverse=False,as_cmap=True)
cm4=cmap=clt.stitch_colormaps(cm2,cm3,name='divergent')
plt.register_cmap(name=cm4.name, cmap=cm4)
#EO.plot_legend('divergent',uni_vmax[0],uni_vmin[0],'test', '%s%swood_gCm2_MOHC.tif' % (path2output,prefix));plt.show()
cm5=clt.cmap_to_perceived_luminosity('divergent')
plt.register_cmap(name=cm5.name, cmap=cm5)
plt.set_cmap('foliage')

"""
#===============================================================================
PART A: DEFINE PATHS AND LOAD IN DATA
#-------------------------------------------------------------------------------
"""
path2data = '/exports/csce/datastore/geos/groups/gcel/CARDAMOM_Brazil/cmip6_scenarios/'
path2output = '/exports/csce/datastore/geos/users/dmilodow/EOlaboratory/EOlab/F2020_CARDAMOM_Brazil/'
prefix = 'Brazil_'
variables = ['wood_gCm2','som_gCm2']
ssps = ['ssp119','ssp126','ssp370','ssp434','ssp585']
uni_vmax = [300,500]
uni_vmin = [0,0]
div_vmax = [2,2]
div_vmin = [-2,-2]
for vv,var in enumerate(variables):
    for ss,ssp in enumerate(ssps):
        var_init = xr.open_rasterio('%s%s%s_initial_MOHC_%s.tif' % (path2data,prefix,var,ssp))[0]
        var_2050 = xr.open_rasterio('%s%s%s_2050_MOHC_%s.tif' % (path2data,prefix,var,ssp))[0]
        var_2050_change = xr.open_rasterio('%s%s%s_2050_change_MOHC_%s.tif' % (path2data,prefix,var,ssp))[0]
        var_2050_change_upper = xr.open_rasterio('%s%s%syr_2050_change_upper_MOHC_%s.tif' % (path2data,prefix,var,ssp))[0]
        var_2050_change_lower = xr.open_rasterio('%s%s%syr_2050_change_lower_MOHC_%s.tif' % (path2data,prefix,var,ssp))[0]
        # apply nodata value
        var_init.values[var_init.values<-3*10**38] = np.nan
        var_2050.values[var_2050.values<-3*10**38] = np.nan
        var_2050_change.values[var_2050_change.values<-3*10**38] = np.nan
        var_2050_change_upper.values[var_2050_change_upper.values<-3*10**38] = np.nan
        var_2050_change_lower.values[var_2050_change_lower.values<-3*10**38] = np.nan

        var_init.values[var_init.values==-9999] = np.nan
        var_2050.values[var_2050.values==-9999] = np.nan
        var_2050_change.values[var_2050_change.values==-9999] = np.nan
        var_2050_change_upper.values[var_2050_change_upper.values==-9999] = np.nan
        var_2050_change_lower.values[var_2050_change_lower.values==-9999] = np.nan

        # generate a mask for significant change at provided confidence level
        mask = np.all((var_2050_change_lower.values<=0,var_2050_change_upper.values>=0),axis=0)

        # convert to Mg / ha
        var_init.values *= 10**4/10**6
        var_2050.values *= 10**4/10**6
        var_2050_change.values *= 10**4/10**6

        # delete existing dataset if present
        prefixes = ['%s%s%s_initial_MOHC_%s.tif' % (path2output,prefix,var,ssp),
                    '%s%s%s_2050_MOHC_%s.tif' % (path2output,prefix,var,ssp),
                    '%s%s%s_2050_change_MOHC_%s' % (path2output,prefix,var,ssp)]
        for pf in prefixes:
            if '%s_data.tif' % pf in os.listdir(path2output):
                os.system("rm %s" % ('%s_data.tif' % (pf)))
            if '%s_display.tif' % pf in os.listdir(path2output):
                os.system("rm %s" % ('%s_display.tif' % (pf)))
            if '%s_confidence_data.tif' % pf in os.listdir(path2output):
                os.system("rm %s" % ('%s_confidence_data.tif' % (pf)))
            if '%s_confidence_display.tif' % pf in os.listdir(path2output):
                os.system("rm %s" % ('%s_confidence_display.tif' % (pf)))

        EO.write_xarray_to_display_layer_GeoTiff(var_init, '%s' % prefixes[0],
                                                'foliage', uni_vmax[vv], uni_vmin[vv])

        EO.write_xarray_to_display_layer_GeoTiff(var_2050, '%s' % prefixes[1],
                                                'foliage', uni_vmax[vv], uni_vmin[vv])

        EO.write_xarray_to_display_layer_GeoTiff(var_2050_change, '%s' % prefixes[2],
                                                'divergent', div_vmax[vv], div_vmin[vv])

        EO.write_xarray_to_display_layer_confidence_levels_GeoTiff(var_2050_change, mask, '%s' % prefixes[2],
                                                'divergent', div_vmax[vv], div_vmin[vv])

EO.plot_legend('foliage',uni_vmax[0],uni_vmin[0],'C$_{wood}$ / Mg C ha$^{-1}$', '%s%swood_gCm2_MOHC.tif' % (path2output,prefix),extend='max')
EO.plot_legend('foliage',uni_vmax[0],uni_vmin[0],'C$_{SOM}$ / Mg C ha$^{-1}$', '%s%ssom_gCm2_MOHC.tif' % (path2output,prefix),extend='max')
EO.plot_legend('divergent',div_vmax[0],div_vmin[0],'C$_{wood}$ / Mg C ha$^{-1}$ yr$^{-1}$', '%s%swood_gCm2_change_MOHC.tif' % (path2output,prefix),extend='both')
EO.plot_legend('divergent',div_vmax[1],div_vmin[1],'C$_{SOM}$ / Mg C ha$^{-1}$ yr$^{-1}$', '%s%ssom_gCm2_change_MOHC.tif' % (path2output,prefix),extend='both')
