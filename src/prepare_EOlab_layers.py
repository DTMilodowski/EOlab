import numpy as np
from matplotlib import pyplot as plt
from osgeo import gdal
import os
import osr
import sys
from netCDF4 import Dataset


import matplotlib as mpl
import matplotlib.cm as cm

# This is a super simple function that loads in a NetCDF file and pulls out the important coordinate
# system info that is needed for writing georeferenced GeoTIFFs.  Since NetCDF files will vary in
# terms of their dimensionality and number of variables, subsequent processing to GeoTIFF layers has
# been separated out into a separate function
def load_NetCDF(NetCDF_file,lat_var = 'lat', lon_var = 'lon'):
    dataset = Dataset(NetCDF_file)
    
    # Get the spatial information from the layer
    Lat = dataset.variables[lat_var][:]
    Long = dataset.variables[lon_var][:]
    
    DataResX = np.abs(Long[0]-Long[1])
    DataResY = np.abs(Lat[0]-Lat[1])
    
    XMinimum = np.min(Long)-DataResX/2.
    YMinimum = np.min(Lat)-DataResY/2.

    geoTransform = [ XMinimum, DataResX, 0, YMinimum, 0, DataResY ]
    
    return dataset, geoTransform

# Convert a python array with float variables into a three-band RGB array with the colours specified
# according to a given colormap and upper and lower limits to the colormap 
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

# Function to write an EO lab data layer from an array
def write_array_to_data_layer_GeoTiff(array,geoTrans, OUTFILE_prefix, EPSG_CODE='4326'):
    NBands = 1
    NRows = 0
    NCols = 0

    # Get array dimensions and flip so that it plots in the correct orientation on GIS platforms
    if len(array.shape) < 2: 
        print 'array has less than two dimensions! Unable to write to raster'
        sys.exit(1)  
    elif len(array.shape) == 2:
        (NRows,NCols) = array.shape
        array = np.flipud(array)
    elif len(array.shape) == 3:
        (NRows,NCols,NBands) = array.shape
        for i in range(0,NBands):
            array[:,:,i] = np.flipud(array[:,:,i])
    else:
        print 'array has too many dimensions! Unable to write to raster'
        sys.exit(1)  
    
        
    # Write GeoTiff
    driver = gdal.GetDriverByName('GTiff')
    driver.Register()

    # set all the relevant geospatial information
    dataset = driver.Create( OUTFILE_prefix+'_data.tif', NCols, NRows, NBands, gdal.GDT_Float32 )
    dataset.SetGeoTransform( geoTrans )
    srs = osr.SpatialReference()
    srs.SetWellKnownGeogCS( 'EPSG:'+EPSG_CODE )
    dataset.SetProjection( srs.ExportToWkt() )
    # write array
    dataset.GetRasterBand(1).SetNoDataValue( -9999 )
    dataset.GetRasterBand(1).WriteArray( array )
    dataset = None
    return 0


# This function is similar to before, except that now it writes two GeoTIFFs - a data layer and a
# display layer.  For the moment, this can only accept a 2D input array
def write_array_to_display_layer_GeoTiff(array, geoTrans, OUTFILE_prefix, cmap, ulim, llim, EPSG_CODE_DATA='4326', ESPG_CODE_DISPLAY='3857'):

    NBands = 1
    NBands_RGB = 3
    NRows = 0
    NCols = 0

    # Get array dimensions and flip so that it plots in the correct orientation on GIS platforms
    if len(array.shape) < 2: 
        print 'array has less than two dimensions! Unable to write to raster'
        sys.exit(1)  
    elif len(array.shape) == 2:
        (NRows,NCols) = array.shape
        array = np.flipud(array)
    else:
        print 'array has too many dimensions! Unable to write to raster'
        sys.exit(1)  

    # Convert RGB array
    rgb_array = convert_array_to_rgb(array,cmap,ulim,llim)
    
    
    # Write Data Layer GeoTiff
    driver = gdal.GetDriverByName('GTiff')
    driver.Register()

    # set all the relevant geospatial information
    dataset = driver.Create( OUTFILE_prefix+'_data.tif', NCols, NRows, NBands, gdal.GDT_Float32 )
    dataset.SetGeoTransform( geoTrans )
    srs = osr.SpatialReference()
    srs.SetWellKnownGeogCS( 'EPSG:'+EPSG_CODE_DATA )
    dataset.SetProjection( srs.ExportToWkt() )
    # write array
    dataset.GetRasterBand(1).SetNoDataValue( -9999 )
    dataset.GetRasterBand(1).WriteArray( array )
    dataset = None

    # Write Display Layer GeoTiff
    driver = gdal.GetDriverByName('GTiff')
    driver.Register()

    # set all the relevant geospatial information
    dataset = driver.Create( 'temp.tif', NCols, NRows, NBands_RGB, gdal.GDT_Float32 )
    dataset.SetGeoTransform( geoTrans )
    srs = osr.SpatialReference()
    srs.SetWellKnownGeogCS( 'EPSG:'+EPSG_CODE_DATA )
    dataset.SetProjection( srs.ExportToWkt() )
    # write array
    dataset.GetRasterBand(1).SetNoDataValue( -9999 )
    dataset.GetRasterBand(1).WriteArray( array )
    dataset = None

    # now use gdalwarp to reproject 
    os.system("gdalwarp -t_srs EPSG:" + EPSG_CODE_DISPLAY + " -srcnodata -9999 -dstnodata -9999 temp.tif " + OUTFILE_prefix+'_display.tif')
    os.system("rm temp.tif")    
    return 0


# A function to produce a simple map legend for quantitative data layers
def plot_legend(cmap,ulim,llim):
    return 0
