# Created By: Susan Meerdink
# 3/27/17
# This python code is the main script for development of spectral libraries
# --------------------------------------------------------------------------------------------------------------------
# Inputs
polyLocation = 'F:\\Classification-Products\\2016_08_15_reference_polygons_withMeta'
metaLocation = 'F:\\Classification-Products\\2016_08_15_reference_polygons_catalog.csv'
dirLocation = 'F:\\Image-To-Image-Registration\\AVIRIS\\*'  # File location for AVIRIS or AVIRIS+MASTER that contains all flightlines
outFigLoc = 'F:\\Dropbox\\Analysis\\Chapter 1 Analysis\\CDALDA_Figures\\'  # Ouput file location for figures generated from analysis
dateTag = '130411'

# Import Modules
import glob
import numpy as np
import rasterio.plot
from rasterio import features
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt
import fiona
from descartes import PolygonPatch
from shapely import *
from shapely.geometry import shape, MultiPolygon
from shapely.ops import transform
from mpl_toolkits.basemap import Basemap
import pyproj
from functools import partial

# Open shapefile and transform to appropriate coordinate system
polyCRS = fiona.open(polyLocation + '.shp', "r")  # Open polygon to get crs information for transformation
polyOriginal = MultiPolygon([shape(pol['geometry']) for pol in fiona.open(polyLocation + '.shp')])  # Open polygon as a MultiPolygon for processing purposes
project = partial(  # Define function for projection process
    pyproj.transform,
    pyproj.Proj(polyCRS.crs),  # source coordinate system
    pyproj.Proj(init='epsg:32611'))  # destination coordinate system, UTM Zone 11 WGS 84
polygons = transform(project, polyOriginal)  # apply projection
print polygons
print len(polygons), "Polygons Found"

# Plot reference polygons
# cm = plt.get_cmap('RdBu')
# num_colours = len(polygons)
# fig = plt.figure()
# ax = fig.add_subplot(111)
# minx, miny, maxx, maxy = polygons.bounds
# w, h = maxx - minx, maxy - miny
# ax.set_xlim(minx - 0.2 * w, maxx + 0.2 * w)
# ax.set_ylim(miny - 0.2 * h, maxy + 0.2 * h)
# ax.set_aspect(1)
# patches = []
# for idx, p in enumerate(polygons):
#     # colour = cm(1. * idx / num_colours)
#     # patches.append(PolygonPatch(p, fc=colour, ec='#555555', alpha=1., zorder=1))
#     patches.append(PolygonPatch(p, fc='#000000', ec='#000000', alpha=1., zorder=1))
# ax.add_collection(PatchCollection(patches, match_original=True))
# ax.set_xticks([])
# ax.set_yticks([])
# plt.title("Plant Species Reference Polygons")
# plt.savefig('F:\\Image-To-Image-Registration\\AVIRIS\\test.png', alpha=True, dpi=300)
# plt.show()
# plt.clf()

# Find image files that spectra will be extracted from
flList = ['FL04']  # 'FL02', 'FL03', 'FL04', 'FL05', 'FL06', 'FL07', 'FL08', 'FL09', 'FL10', 'FL11'
imageList = []
for fl in flList:  # loop through flightline files to find specific date
    imageLocation = dirLocation + fl + '\\6 - Spectral Correction Files\\*' + dateTag + '*'
    for name in glob.glob(imageLocation):  # Ignore header files
        if '.hdr' not in name:
            print 'Found', name
            imageList.append(name)  # add file to list of image files for future processing

spectralLibrary = []
for i in imageList:  # loop through flightline files to extract spectra
    print 'Extracting spectra from', i
    imgFile = rasterio.open(i, 'r')  # Open raster image
    iName = i.split('\\')[-1]  # Get file name

    # Plot reference polygons on top of current image
    fig, axMap = plt.subplots(num=None, figsize=(4, 3), dpi=300, facecolor='w', edgecolor='k')
    rasterio.plot.show((imgFile, 137), ax=axMap, cmap='hot', clim=(1, 10000))  # Show band 137 of current image
    patches = []
    for idx, p in enumerate(polygons):
        patches.append(PolygonPatch(p, fc='#FFFFFF', ec='#FFFFFF', alpha=1., zorder=1))
    axMap.add_collection(PatchCollection(patches, match_original=True))
    plt.xticks(fontsize=6)
    plt.yticks(fontsize=6)
    plt.title("Plant Species Reference Polygons with \n" + iName, fontsize=6)
    plt.savefig(outFigLoc + iName + '_reference_poly.png', alpha=True, dpi=300)
    plt.clf()

    # Extract spectra from image
    # Note: cannot read entire AVIRIS file as it is too large, need to loop by band
    for bandNum in range(1, 2):  # Loop through bands (note: bands are indexed from 1)
        imgData = imgFile.read(bandNum)

        for idx in range(0, 1):  # Loop through polygons len(polygons)
            print polygons[idx]['geometry']
            print polygons[idx]['properties']['Polygon_ID']
            # Create Mask that has 1 for locations with polygons and 0 for non polygon locations
            polygonMask = rasterio.features.geometry_mask(polyCRS[idx]['geometry'], out_shape=imgFile.shape, transform=imgFile.transform, all_touched=True)
            fig, axMap = plt.subplots(num=None, figsize=(4, 3), dpi=300, facecolor='w', edgecolor='k')
            rasterio.plot.show(polygonMask, ax=axMap)

            # create a masked numpy array
            masked_array = np.ma.array(data=imgData, mask=polygonMask.astype(bool))
            print (masked_array > 0).nonzero()
            masked_data = masked_array[(masked_array > 0).nonzero()]
            print masked_data

