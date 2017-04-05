# Created By: Susan Meerdink
# 3/27/17
# This python code is the main script for development of spectral libraries
# --------------------------------------------------------------------------------------------------------------------
# Inputs
polyLocation = 'F:\\Classification-Products\\2016_08_15_reference_polygons_withMeta.shp'
metaLocation = 'F:\\Classification-Products\\2016_08_15_reference_polygons_catalog.csv'
dirLocation = 'F:\\Image-To-Image-Registration\\AVIRIS\\*' #File location for AVIRIS or AVIRIS+MASTER that contains all flightlines
dateTag = '130411'

# Import Modules
import glob
import numpy as np
import rasterio
import ogr
import matplotlib.pyplot as plt
import fiona

# Open shapefile
with fiona.open(polyLocation, "r") as shapefile:
    polygon = [feature["geometry"] for feature in shapefile]

# Find image files that spectra will be extracted from
flList = ['FL02', 'FL03', 'FL04', 'FL05', 'FL06', 'FL07', 'FL08', 'FL09', 'FL10', 'FL11']
imageList = []
for fl in flList: #loop through flightline files to find specific date
    imageLocation = dirLocation + fl + '\\6 - Spectral Correction Files\\*' + dateTag + '*'
    for name in glob.glob(imageLocation): #Ignore header files
        if '.hdr' not in name:
            print 'Found', name
            imageList.append(name) #add file to list of image files for future processing

for i in imageList: #loop through flightline files to extract spectra
    img = rasterio.open(i)
    print img.shape

    # Extract spectra from image

