"""
PotentialAGB_app.py
This function creates the required layers for potential AGB EO lab applications.
Specific layers produced:
    - Observed AGB
    - Potential AGB
    - Deficit
    - Sequestration Potential
    - Grid cell areas (for on-the-fly observations)
Both data (WGS84) and display (RGB, Web Mercator) layers are produced, with an
associated legend for the latter.

Input layers are assumed to be netcdf files produced with the scripts in the
PotentialBiomassRFR library.

Input arguments help to point to the correct files for the country of interest.
If the "country" option is specified as "None", then by default no clipping of
the data will be undertaken, otherwise the data will be clipped to the national
boundaries of the specified country.
"""

## INPUT ARGUMENTS
country = 'Colombia'
DATADIR = '/home/dmilodow/DataStore_DTM/FOREST2020/PotentialBiomassRFR/output/'
SAVEDIR = '/home/dmilodow/DataStore_DTM/EOlaboratory/EOlab/ColombiaPotentialAGB/'
NetCDF_file = 'COL_003_AGB_potential_RFR_worldclim_soilgrids_2pass.nc'


## IMPORT PACKAGES
from import_pkgs import *

## READ IN DATA

sf = shapefile.Reader("/home/dmilodow/DataStore_DTM/EOlaboratory/Areas/NaturalEarth/10m_cultural/ne_10m_admin_0_countries")

vertices = []
for shape_rec in sf.shapeRecords():
    if shape_rec.record[3] == country_name:
        pts = shape_rec.shape.points
        prt = list(shape_rec.shape.parts) + [len(pts)]
        for i in range(len(prt) - 1):
            vertices.append([])
            for j in range(prt[i], prt[i+1]):
                vertices[i].append((pts[j][0], pts[j][1]))


plt.figure(1, facecolor='White',figsize=[2, 1])
plt.show()

DATADIR = '/home/dmilodow/DataStore_DTM/FOREST2020/PotentialBiomassV2/output/'
SAVEDIR = '/home/dmilodow/DataStore_DTM/EOlaboratory/EOlab/ColombiaPotentialAGB/'
NetCDF_file = 'COL_v002_AGBpot_mean_WC2_SOILGRIDS_GridSearch.nc'

ds,geoTrans = EO.load_NetCDF(DATADIR+NetCDF_file,lat_var = 'lat', lon_var = 'lon')
lat = np.asarray(ds.variables['lat'])
lon = np.asarray(ds.variables['lon'])
vars = ['AGB_mean','AGBpot_mean','training','areas']

# create clipping mask
S=-5.
N=13.
E= -66.5
W=-79.05
lat_mask = np.all((lat<=N,lat>=S),axis=0)
lon_mask = np.all((lon<=E,lon>=W),axis=0)
n_lat = lat_mask.sum()
n_lon = lon_mask.sum()
lat_clip = lat[lat_mask]
lon_clip = lon[lon_mask]
array_mask = np.ix_(lat_mask,lon_mask)

latgrid = np.zeros((n_lat,n_lon))
longrid = np.zeros((n_lat,n_lon))
for ilat in range(0,n_lat):
  longrid[ilat,:] = lon_clip.copy()

for ilon in range(0,n_lon):
  latgrid[:,ilon] = lat_clip.copy()

# now loop through all polygons, and find the pixels which fall within them
country_mask = np.zeros((n_lat,n_lon))

for vv in range(0,len(vertices)):
      temp1,temp2,inside = aux.points_in_poly(longrid.ravel(),latgrid.ravel(),vertices[vv])
      country_mask[inside.reshape(n_lat,n_lon)] = 1

# apply clip
dataset = {}
for vv in range(0,len(vars)):
    dataset[vars[vv]]= np.asarray(ds.variables[vars[vv]])[array_mask]
    dataset[vars[vv]][country_mask==0]=-9999


# create new geoTrans object to account for clipping
geoTrans[0] = np.min(lon_clip)
if geoTrans[5]>0:
  geoTrans[3]=np.min(lat_clip)
else:
  geoTrans[3]=np.max(lat_clip)

# sequestration potential is defined by pixels with positive potential biomass that
# are not already forests
dataset['seqpot_mean'] = dataset['AGBpot_mean']-dataset['AGB_mean']
#dataset['seqpot_mean'][dataset['forests']==1] = 0.
dataset['seqpot_mean'][dataset['seqpot_mean']<0] = 0.
dataset['seqpot_mean'][dataset['AGB_mean']==-9999] = -9999.

dataset['training'][dataset['training']<1] = -9999.

vars = ['AGB_mean','AGBpot_mean','seqpot_mean','training']
cmaps = ['viridis','viridis','plasma','viridis']
ulims = [400.,400.,200.,1.]
llims = [0.,0.,0.,0.]
axis_labels = ['AGB$_{obs}$ / Mg(C) ha$^{-1}$', 'AGB$_{potential}$ / Mg(C) ha$^{-1}$', 'Sequestration potential / Mg(C) ha$^{-1}$', 'Forest mask (1 = Forest)']

for vv in range(0,len(vars)):
    print(vars[vv])
    file_prefix = SAVEDIR + 'colombia_' + vars[vv]

    # delete existing dataset if present
    if 'colombia_'+vars[vv]+'_data.tif' in os.listdir(SAVEDIR):
        os.system("rm %s" % (SAVEDIR+'colombia_'+vars[vv]+'_data.tif'))
    if 'colombia_'+vars[vv]+'_display.tif' in os.listdir(SAVEDIR):
        os.system("rm %s" % (SAVEDIR+'colombia_'+vars[vv]+'_display.tif'))

    EO.write_array_to_display_layer_GeoTiff(dataset[vars[vv]], geoTrans, file_prefix, cmaps[vv], ulims[vv], llims[vv])
    if vars[vv]!='training':
        EO.plot_legend(cmaps[vv],ulims[vv],llims[vv],axis_labels[vv], file_prefix)

"""
rows, cols = dataset[vars[0]].shape
latitude = np.arange(geoTrans[3],rows*geoTrans[5]+geoTrans[3],geoTrans[5])[:rows]
longitude =  np.arange(geoTrans[0],cols*geoTrans_rs[1]+geoTrans_rs[0],geoTrans_rs[1])[:cols]
areas = geo.calculate_cell_area_array(latitude,longitude, area_scalar = 1./10.**4,cell_centred=False)
"""

# loop through the variables, multiplying by cell areas to give values in Mg
for vv in range(0,len(vars)):
    print(vars[vv])

    if 'colombia_'+vars[vv]+'_total_data.tif' in os.listdir(SAVEDIR):
        os.system("rm %s" % (SAVEDIR+'colombia_'+vars[vv]+'_total_data.tif'))

    file_prefix = SAVEDIR + 'colombia_' + vars[vv] + '_total'

    out_array = dataset[vars[vv]] * dataset['areas']
    out_array[dataset[vars[vv]]==-9999]=-9999
    EO.write_array_to_data_layer_GeoTiff(out_array, geoTrans, file_prefix)
    out_array=None

# Also want to write cell areas to file.  However, as this will be compared against other layers, need to carry across
# nodata values
areas_out = dataset['areas']
areas_out[np.asarray(dataset[vars[0]])==-9999] = -9999
if 'colombia_cell_areas_data.tif' in os.listdir(SAVEDIR):
    os.system("rm %s" % (SAVEDIR+'colombia_cell_areas_data.tif'))
area_file_prefix = SAVEDIR + 'colombia_cell_areas'
EO.write_array_to_data_layer_GeoTiff(areas_out, geoTrans, area_file_prefix)
