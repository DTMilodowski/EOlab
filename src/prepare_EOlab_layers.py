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


def write_array_to_data_layer_GeoTiff(array,<details>, OUTFILE, EPSG_CODE):

    # first set all the relevant geospatial information
    dataset_ds = driver.Create( 'temp.tif', NCols, NRows, NBands, gdal.GDT_Float32 )
    dataset_ds.SetGeoTransform( [ XMinimum, DataResX, 0, YMinimum, 0, DataResY ] )
    srs = osr.SpatialReference()
    srs.SetWellKnownGeogCS( 'EPSG:'+EPSG_CODE_INIT )
    dataset_ds.SetProjection( srs.ExportToWkt() )
    # write array
    dataset_ds.GetRasterBand(1).SetNoDataValue( -9999 )
    dataset_ds.GetRasterBand(1).WriteArray( Ctot_u )
    dataset_ds = None
    # now use gdalwarp to reproject 
    os.system("gdalwarp -t_srs EPSG:" + EPSG_CODE_TARGET + " -srcnodata -9999 -dstnodata -9999 temp.tif " + OUTFILE)
    os.system("rm temp.tif")    

    return 0

def write_array_to_display_layer_GeoTiff(array,<details>):
    return 0
