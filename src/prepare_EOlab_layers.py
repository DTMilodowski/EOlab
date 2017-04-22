import numpy as np
from matplotlib import pyplot as plt
from osgeo import gdal
import os
import osr
import sys
from netCDF4 import Dataset


import matplotlib as mpl
import matplotlib.cm as cm

def convert_array_to_rgb(array, cmap, ulim, llim):
  norm = mpl.colors.Normalize(vmin=llim, vmax=ulim)
  
  rows, cols = array.shape
  rgb_array = np.zeros((rows,cols,3),dtype=np.uint8)*np.nan
  for i in range(0,rows):
    for j in range(0,cols):
      if np.isfinite(array[i,j]):
        rgb_array[i,j,:]= cm.ScalarMappable(norm=norm,cmap=cmap).to_rgba(array[i,j])
      else:
        rgb_array[i,j,:] = np.asarray([255,0,255]) # check this
  
  return rgb_array


