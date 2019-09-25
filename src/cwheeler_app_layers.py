# Import general libraries
import numpy as np
import xarray as xr
# import plotting libraries
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()
# import custom libraries
import prepare_EOlab_layers as EO
# define cmaps that might be used (examples from seaborn library)
plt.register_cmap(name='divergent', cmap=sns.diverging_palette(275,150,l=66,s=90,as_cmap=True))
plt.register_cmap(name='seagreen', cmap=sns.light_palette('seagreen',as_cmap=True))
"""
#===============================================================================
PART A: DEFINE PATHS AND LOAD IN DATA; SET UP PRESENTATION SPECS
Change according to your directory structure and other requirements
#-------------------------------------------------------------------------------
"""
raster_file = 'JAL_AGB_2018_WGS84_4362.tif'
path2data  = '/home/dmilodow/DataStore_DTM/FOREST2020/CWheeler/'
path2output = '/home/dmilodow/DataStore_DTM/FOREST2020/CWheeler/'

# load potential biomass models from netdf file
dataset = xr.open_rasterio('%s%s' % (path2data, raster_file))[0]
dataset.values[dataset.values>3e38] = np.nan # convert nodata value to -9999. This might change for other datasets

# set some variables for presentation e.g. colormap specs
file_prefix = '%s%s' % (path2output,raster_file.split('.')[0]) # a prefix for output files (full path included)
cmap = 'seagreen' # change as required
axis_label = 'AGB / Mg ha$^{-1}$' # change as required
ulim = 150 # upper limit of color ramp
llim = 0 # lower limit of color ramp

# check that things look ok (this will just plot to screen).
# Note that aspect ratio might not be preserved, but this isn't important
# Main point is to check that upper and lower limits of the colormap are
# appropriate, and that you are happy with the colormap
dataset.plot(cmap=cmap,vmin=llim,vmax=ulim)
plt.show()
"""
#===============================================================================
PART B: CREATE DISPLAY LAYER
#-------------------------------------------------------------------------------
"""
EO.write_xarray_to_display_layer_GeoTiff(dataset, file_prefix, cmap, ulim, llim) # this will create both a display and a data layer
EO.plot_legend(cmap,ulim,llim,axis_label, file_prefix) # this will create a png file with the legend, that displays nicely on the app
