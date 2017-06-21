
# coding: utf-8

# In[116]:

import pandas as pd
import numpy as np
import xarray as xr
import datacube
from datacube.helpers import ga_pq_fuser
from datacube.api import make_mask
from affine import Affine
import fiona
import rasterio
import rasterio.features
import matplotlib.pyplot as plt


# In[117]:

# Parameters for the extract
startDate = '1995-01-01'
numMonths = 180
shapefileName = '/home/147/ss8137/Documents/koala_obs2/tile222.shp' 


# In[118]:

def load_multiple_masked_dataset(platforms, query, bands):
    nbars = []
    fcs = []
    for platform in platforms:
        product_name = '{}_{}_albers'.format(platform, 'nbar')
        ds = dc.load(product=product_name, measurements=bands, group_by='solar_day', **query)
        
        if not ds:
            continue
        product_fc= '{}_{}_albers'.format(platform, 'fc')
        fcov = dc.load(product = product_fc, measurements=pfc_measurements, group_by='solar_day', **query)
        
        if not fcov:
            continue
        mask_product = '{}_{}_albers'.format(platform, 'pq')
        sensor_pq = dc.load(product=mask_product, fuse_func=ga_pq_fuser, group_by='solar_day', **query)
        if not sensor_pq:
            continue
        cloud_free = make_mask(sensor_pq.pixelquality, ga_good_pixel=True)
        ds = ds.where(cloud_free)
        ds['product'] = ('time', np.repeat(product_name, ds.time.size))
        fcov = fcov.where(cloud_free)
        fcov['product'] = ('time', np.repeat(product_fc, fcov.time.size))
        nbars.append(ds)
        fcs.append(fcov)

    if len(nbars)>0:
        combined_nbar = xr.concat(nbars, dim='time')
        combined_nbar = combined_nbar.isel(time=combined_nbar.time.argsort())
        return combined_nbar
    elif len(fcov)>0:
        combined_fc = xr.concat(fcs, dim='time')
        combined_fc = combined_fc.isel(time=combined_fc.time.argsort())
        return combined_fc
    else:
        return None


# In[119]:

dc = datacube.Datacube(app='point-ndvi-savi-recipe')


# In[120]:

# Generate a list of monthly dates
dateList = pd.date_range(startDate, periods=numMonths+1, freq='1M').tolist()
print("Processing data from " + str(dateList[0]) + " to " + str(dateList[-1]))


# In[121]:

# Read in the points. This has to be EPSG:3577. 
with fiona.open(shapefileName, "r") as shapefile:
    # Get the total extent of all the data in the shapefile
    bounds = shapefile.bounds 
    # Extract the individual geometries
    geoms = [feature["geometry"] for feature in shapefile]
    numGeoms = len(geoms)
    # You could extract anything here - I'm just grabbing a list of the ID fields
    ids = [feature["id"] for feature in shapefile]
    print("Order of Geometries is" + str(ids))


# In[122]:

# Initialise Arrays

NDVI = np.zeros((numMonths,numGeoms))
SAVI = np.zeros((numMonths,numGeoms))


# In[123]:

# Use 5-7 Landsats.
platforms = ['ls5', 'ls7']
# Only Want the bands to calculate NDVI
bands_of_interest = ['nir','red']

pfc_measurements=['BS']


# In[124]:

# Loop through the months. 
for i in range(numMonths):
    print(dateList[i])
    # Setup The Query using the bounds of the shapefile and a 1 month interval
    query = {'time': (dateList[i],dateList[i+1]),'x': (bounds[0], bounds[2]),'y': (bounds[1], bounds[3]) ,'crs':'EPSG:3577'}
    # Load all the data for the month from the DataCube
    combined_nbar = load_multiple_masked_dataset(platforms, query, bands_of_interest)
    if combined_nbar is not None:
        # Compute the per pixel NDVI  and SAVI over the period
        #lfac = combined_fc.BS
        monthlyNDVI = ((combined_nbar.nir - combined_nbar.red) / (combined_nbar.nir + combined_nbar.red))
        # Extract each of the geometries 
        for j in range(numGeoms):
            # Build an analysis mask using a crazy 1-liner. Doing this every loop is dumb but it's the premature optimization + redability thing. Hey, it works
            mask=rasterio.features.geometry_mask([geoms[j]],out_shape=combined_nbar.geobox.shape,transform=combined_nbar.geobox.affine,all_touched=False,invert=True)
            # Extract the NDVI
            NDVI[i,j] = monthlyNDVI.where(mask)
            # Count the number of pixels of ndvi in the mask area
            countOfPixels_ndvi = np.sum(np.isfinite(monthlyNDVI.where(mask)))
            
    combined_fc = load_multiple_masked_dataset(platforms, query, bands_of_interest)
    if combined_fc is not None:
        lfac = combined_fc.BS
        monthlySAVI = ((combined_nbar.nir - combined_nbar.red) / (combined_nbar.nir + combined_nbar.red + lfac))*(1 + lfac)
        for j in range(numGeoms):
        # Build an analysis mask using a crazy 1-liner. Doing this every loop is dumb but it's the premature optimization + redability thing. Hey, it works
            mask=rasterio.features.geometry_mask([geoms[j]],out_shape=combined_fc.geobox.shape,transform=combined_fc.geobox.affine,all_touched=False,invert=True)
            # Extract the SAVI
            SAVI[i, j] = monthlySAVI.where(mask)
            # Count the number of pixels of savi in the mask area
            countOfPixels_savi = np.sum(np.isfinite(monthlySAVI.where(mask)))
    else:
        # There is noData. 
        SAVI[i,:] = np.nan
#here is the error
ValueError                                Traceback (most recent call last)
<ipython-input-133-6abff40978eb> in <module>()
     16             mask=rasterio.features.geometry_mask([geoms[j]],out_shape=combined_nbar.geobox.shape,transform=combined_nbar.geobox.affine,all_touched=False,invert=True)
     17             # Extract the NDVI
---> 18             NDVI[i,j] = monthlyNDVI.where(mask)
     19             # Count the number of pixels of ndvi in the mask area
     20             countOfPixels_ndvi = np.sum(np.isfinite(monthlyNDVI.where(mask)))

ValueError: setting an array element with a sequence.

# In[ ]:

# Save the results to a file
np.savetxt("/home/147/ss8137/Documents/koala_obs2/NDVI_t222.csv", NDVI, delimiter=',')
np.savetxt("/home/147/ss8137/Documents/koala_obs2/SAVI_t222.csv", SAVI, delimiter=',')


# In[ ]:

# Plot the results in a plot. Actually, 2 plots. Well, 2 subplots in 1 plot. Probaby just one plot I guess if you want to be specific.
# This would be cooler with the xarray magic plottingness if I could work out how to make  NDVI etc as xarrays. 
f, (ax1, ax2) = plt.subplots(2, 1,figsize=(16,8),sharex=True)
# Plot the two metrics
ax1.plot(dateList[1:],NDVI)
ax2.plot(dateList[1:],SAVI)
# Add grids because grids look scientific
ax1.grid(True)
ax2.grid(True)
# Add in titles and labels after the reviewer complained
ax1.set_title(' NDVI')
ax2.set_title('SAVI')
ax2.set_xlabel('Date')
ax1.set_ylabel('NDVI')
ax2.set_ylabel('SAVI')
# Add a legend. Or two.
ax1.legend(ids, loc='upper left')
ax2.legend(ids, loc='upper left')


# In[ ]:



