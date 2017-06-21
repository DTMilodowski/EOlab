import numpy as np
from matplotlib import pyplot as plt
from osgeo import gdal
import os
import osr
import sys
from netCDF4 import Dataset


import matplotlib as mpl
import matplotlib.cm as cm
from matplotlib import ticker

import scipy
from scipy import ndimage, signal


from matplotlib import rcParams
# Set up some basiic parameters for the plots
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['arial']
rcParams['font.size'] = 9
rcParams['legend.numpoints'] = 1
axis_size = rcParams['font.size']

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
    YMaximum = np.max(Lat)+DataResY/2.

    #geoTransform = [ XMinimum, DataResX, 0, YMinimum, 0, DataResY ]
    geoTransform = [ XMinimum, DataResX, 0, YMaximum, 0, -DataResY ]
    
    return dataset, geoTransform

# A function to resample an array to a higher resolution.  The resampling is specified using an scalar,
# which should be an integer, and represents the number of divisions each cell is to be split into.
# This only resamples to higher resolutions, and does not deal with aggregating to lower resolutions.
# The main reason for using this is to increase the accuracy of area based queries using polygon
# shapefiles in EO lab applications. vars is a list of variables which you'd like to resample
def resample_dataset(dataset,geoTransform,vars,resampling_scalar):
    ds = {}
    for vv in range(0,len(vars)):
        print vars[vv]
        ds[vars[vv]], geoTrans = resample_array(np.asarray(dataset.variables[vars[vv]]),geoTransform,resampling_scalar)
    return ds, geoTrans

def resample_array(array,geoTransform,resampling_scalar):
    rs = resampling_scalar    
    rows,cols = array.shape
    array_temp = np.zeros((rows*rs,cols*rs))

    # Fill the new array with the original values
    array_temp[::rs,::rs] = array

    # Define the convolution kernel           
    kernel_1d = scipy.signal.boxcar(rs)
    kernel_2d = np.outer(kernel_1d, kernel_1d)

    # Apply the kernel by convolution, seperately in each axis
    array_out = scipy.signal.convolve(array_temp, kernel_2d, mode="same")
    
    #for ii in range(0,rows):
    #    for jj in range(0,cols):
    #        array_out[ii*rs:(ii+1)*rs,jj*rs:(jj+1)*rs]=array[ii,jj]
    geoTrans = [geoTransform[0], geoTransform[1]/float(rs), geoTransform[2], geoTransform[3], geoTransform[4], geoTransform[5]/float(rs)]
    return array_out, geoTrans

# Function to load a GeoTIFF band plus georeferencing information.  Only loads one band,
# which is band 1 by default
def load_GeoTIFF_band_and_georeferencing(File,band_number=1):
    
    driver = gdal.GetDriverByName('GTiff')
    driver.Register()

    try:
        ds = gdal.Open(File)
    except RuntimeError, e:
        print 'unable to open ' + File
        print e
        sys.exit(1)
        
    source_band = ds.GetRasterBand(band_number)
    if source_band is None:
        print "BAND MISSING"
        sys.exit(1)  

    array = np.array(ds.GetRasterBand(band_number).ReadAsArray(),dtype=np.float64)
    geoTrans = ds.GetGeoTransform()
    coord_sys = ds.GetProjectionRef()

    return array, geoTrans, coord_sys


# Convert a python array with float variables into a three-band RGB array with the colours specified
# according to a given colormap and upper and lower limits to the colormap 
def convert_array_to_rgb(array, cmap, ulim, llim, nodatavalue=-9999):
    norm = mpl.colors.Normalize(vmin=llim, vmax=ulim)
    rgb_array= cm.ScalarMappable(norm=norm,cmap=cmap).to_rgba(array)[:,:,:-1]*255
    mask = np.any((~np.isfinite(array),array==nodatavalue),axis=0)

    rgb_array[mask,0]=255.
    rgb_array[mask,1]=0.
    rgb_array[mask,2]=255.
    """  
    rows, cols = array.shape
    rgb_array = np.zeros((rows,cols,3),dtype=np.uint8)*np.nan
    for i in range(0,rows):
        for j in range(0,cols):
            if np.isfinite(array[i,j]):
                rgb_array[i,j,:]= cm.ScalarMappable(norm=norm,cmap=cmap).to_rgba(array[i,j])[:3]
                rgb_array[i,j,:]*=255
            else:
                rgb_array[i,j,:] = np.asarray([255,0,255]) # check this
    """
    return rgb_array

# Function to write an EO lab data layer from an array
def write_array_to_data_layer_GeoTiff(array,geoTrans, OUTFILE_prefix, EPSG_CODE='4326', north_up=True):
    NBands = 1
    NRows = 0
    NCols = 0

    if north_up:
        # for north_up array, need the n-s resolution (element 5) to be negative
        if geoTrans[5]>0:
            geoTrans[5]*=-1
            geoTrans[3] = geoTrans[3]-(array.shape[0]+1.)*geoTrans[5]
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

    else:
        # for north_up array, need the n-s resolution (element 5) to be positive
        if geoTrans[5]<0:
            geoTrans[5]*=-1
            geoTrans[3] = geoTrans[3]-(array.shape[0]+1.)*geoTrans[5]
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
def write_array_to_display_layer_GeoTiff(array, geoTrans, OUTFILE_prefix, cmap, ulim, llim, EPSG_CODE_DATA='4326', EPSG_CODE_DISPLAY='3857', north_up=True):

    NBands = 1
    NBands_RGB = 3
    NRows = 0
    NCols = 0

    if north_up:
        # for north_up array, need the n-s resolution (element 5) to be negative
        if geoTrans[5]>0:
            geoTrans[5]*=-1
            geoTrans[3] = geoTrans[3]-(array.shape[0]+1.)*geoTrans[5]
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

    else:
        # for north_up array, need the n-s resolution (element 5) to be positive
        if geoTrans[5]<0:
            geoTrans[5]*=-1
            geoTrans[3] = geoTrans[3]-(array.shape[0]+1.)*geoTrans[5]
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
    dataset = driver.Create( 'temp.tif', NCols, NRows, NBands_RGB, gdal.GDT_Byte )
    dataset.SetGeoTransform( geoTrans )
    srs = osr.SpatialReference()
    srs.SetWellKnownGeogCS( 'EPSG:'+EPSG_CODE_DATA )
    dataset.SetProjection( srs.ExportToWkt() )
    # write array
    for bb in range(0,3):
        #dataset.GetRasterBand(bb+1).SetNoDataValue( -9999 )
        dataset.GetRasterBand(bb+1).WriteArray( rgb_array[:,:,bb] )
    dataset = None

    # now use gdalwarp to reproject 
    #os.system("gdalwarp -t_srs EPSG:" + EPSG_CODE_DISPLAY + " -srcnodata -9999 -dstnodata -9999 temp.tif " + OUTFILE_prefix+'_display.tif')
    os.system("gdalwarp -t_srs EPSG:" + EPSG_CODE_DISPLAY + " temp.tif " + OUTFILE_prefix+'_display.tif')
    os.system("rm temp.tif")    
    return 0


# A function to produce a simple map legend for quantitative data layers
def plot_legend(cmap,ulim,llim,axis_label, OUTFILE_prefix):
    norm = mpl.colors.Normalize(vmin=llim, vmax=ulim)
    plt.figure(1, facecolor='White',figsize=[2, 1])
    ax = plt.subplot2grid((1,1),(0,0))
    cb = mpl.colorbar.ColorbarBase(ax,cmap=cmap,norm=norm,orientation='horizontal')
    tick_locator = ticker.MaxNLocator(nbins=5)
    cb.locator = tick_locator
    cb.update_ticks()
    cb.set_label(axis_label,fontsize = axis_size)
    plt.tight_layout()
    plt.savefig(OUTFILE_prefix+'_legend.png')
    #plt.show()
    return 0
