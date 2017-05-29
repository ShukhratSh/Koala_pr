
# coding: utf-8

# In[ ]:

get_ipython().magic('pylab notebook')
from __future__ import print_function
import datacube
import xarray as xr
from datacube.storage import masking
from datacube.storage.masking import mask_to_dict
from datacube.helpers import ga_pq_fuser
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot as plt
import matplotlib.dates
import fiona
import shapely
import shapely.geometry
from shapely.geometry import shape
import rasterio


# In[ ]:

dc = datacube.Datacube(app='NDVI,SAVI calculation based on the observed points')


# In[65]:

#This defines the function that converts a linear vector file into a string of x,y coordinates
def geom_query(geom, geom_crs='EPSG:4326'):
    """
    Create datacube query snippet for geometry
    """
    return {
        'x': (geom.bounds[0], geom.bounds[2]),
        'y': (geom.bounds[1], geom.bounds[3]),
        'crs': geom_crs
    }


def warp_geometry(geom, crs_crs, dst_crs):
    """
    warp geometry from crs_crs to dst_crs
    """
    return shapely.geometry.shape(rasterio.warp.transform_geom(crs_crs, dst_crs, shapely.geometry.mapping(geom)))


def transect(data, geom, resolution, method='nearest', tolerance=None):
    """
    
    """
    dist = [i for i in range(0, int(geom.length), resolution)]
    points = list(zip(*[geom.interpolate(d).coords[0] for d in dist]))
    indexers = {
        data.crs.dimensions[0]: list(points[1]),
        data.crs.dimensions[1]: list(points[0])        
    }
    return data.sel_points(xr.DataArray(dist, name='distance', dims=['distance']),
                           method=method,
                           tolerance=tolerance,
                           **indexers)


# In[ ]:

#### DEFINE SPATIOTEMPORAL RANGE AND BANDS OF INTEREST
#Select polyline, replacing /g/... with /*your_path_here*/your_file.shp

vec_fname = '/home/147/ss8137/Documents/koala_obs/koala_obs_1995_2009.shp' 

with fiona.open(vec_fname) as src:
    geom = shape(src[0]['geometry'])

#Define temporal range
start_of_epoch = '1995-01-01'
#need a variable here that defines a rolling 'latest observation'
end_of_epoch =  '2009-12-31'

#Define wavelengths/bands of interest, remove this kwarg to retrieve all bands
bands_of_interest = [#'blue',
                     #'green',
                     'red', 
                     'nir',
                     #'swir1', 
                     #'swir2'
                     ]

#Define sensors of interest
sensors  = ['ls7','ls5']

query = {'time': (start_of_epoch, end_of_epoch),}
query.update(geom_query(geom)) 
query['crs'] = 'EPSG:4326'


# In[ ]:

print (query)


# In[ ]:

#Group PQ by solar day to avoid idiosyncracies of N/S overlap differences in PQ algorithm performance
pq_albers_product = dc.index.products.get_by_name(sensors[0]+'_pq_albers')
valid_bit = pq_albers_product.measurements['pixelquality']['flags_definition']['contiguous']['bits']


# In[ ]:

#Define which pixel quality artefacts you want removed from the results
mask_components = {'cloud_acca':'no_cloud',
'cloud_shadow_acca' :'no_cloud_shadow',
'cloud_shadow_fmask' : 'no_cloud_shadow',
'cloud_fmask' :'no_cloud',
'blue_saturated' : False,
'green_saturated' : False,
'red_saturated' : False,
'nir_saturated' : False,
'swir1_saturated' : False,
'swir2_saturated' : False,
'contiguous':True}


# In[ ]:

#retrieve the NBAR and PQ for the spatiotemporal range of interest

#Retrieve the NBAR and PQ data for sensor n
sensor_clean = {}
for sensor in sensors:
    #Load the NBAR and corresponding PQ
    sensor_nbar = dc.load(product= sensor+'_nbar_albers', group_by='solar_day', measurements = bands_of_interest,  **query)
    sensor_pq = dc.load(product= sensor+'_pq_albers', group_by='solar_day', fuse_func=ga_pq_fuser, **query)
    #grab the projection info before masking/sorting
    crs = sensor_nbar.crs
    crswkt = sensor_nbar.crs.wkt
    affine = sensor_nbar.affine
    #This line is to make sure there's PQ to go with the NBAR
    sensor_nbar = sensor_nbar.sel(time = sensor_pq.time)
    #Apply the PQ masks to the NBAR
    cloud_free = masking.make_mask(sensor_pq, **mask_components)
    good_data = cloud_free.pixelquality.loc[start_of_epoch:end_of_epoch]
    sensor_nbar = sensor_nbar.where(good_data)
    sensor_clean[sensor] = sensor_nbar


# In[66]:


#Concatenate the data from different sensors together and sort by time
nbar_clean = xr.concat(sensor_clean.values(), dim='time')
time_sorted = nbar_clean.time.argsort()
nbar_clean = nbar_clean.isel(time=time_sorted)
nbar_clean.attrs['crs'] = crs
nbar_clean.attrs['affine'] = affine

#Extract the observation data volume
geom_o = warp_geometry(geom, query['crs'], crs.wkt)
obs = transect(nbar_clean, geom_o, 25)


# In[52]:

print('The number of time slices at this location is '+ str(nbar_clean.red.shape[0]))


# In[54]:

#select time slice of interest - this is trial and error until you get a decent image
time_slice_i = 350
rgb = nbar_clean.isel(time =time_slice_i).to_array(dim='color').sel(color=['red', 'nir']).transpose('y', 'x', 'color')
#rgb = nbar_clean.isel(time =time_slice).to_array(dim='color').sel(color=['swir1', 'nir', 'green']).transpose('y', 'x', 'color')
fake_saturation = 4500
clipped_visible = rgb.where(rgb<fake_saturation).fillna(fake_saturation)
max_val = clipped_visible.max(['y', 'x'])
scaled = (clipped_visible / max_val)


# In[56]:

#View the polyline vector on the imagery
fig = plt.figure(figsize =(12,6))
plt.scatter(x=obs.coords['x'], y=obs.coords['y'], c='r') #turn this on or off to show location of transect
plt.imshow(scaled, interpolation = 'nearest',
           extent=[scaled.coords['x'].min(), scaled.coords['x'].max(), 
                   scaled.coords['y'].min(), scaled.coords['y'].max()])

date_ = nbar_clean.time[time_slice_i]
plt.title(date_.astype('datetime64[D]'))
plt.show()


# In[ ]:



