# Created By: Susan Meerdink
# 3/27/17
# This python code is the main script for development of spectral libraries
# --------------------------------------------------------------------------------------------------------------------
# Inputs
polyLocation = 'F:\\Classification-Products\\2016_08_15_reference_polygons_withMeta'
metaLocation = 'F:\\Classification-Products\\2016_08_15_reference_polygons_catalog.csv'
dirLocation = 'F:\\Image-To-Image-Registration\\AVIRIS\\*'  # File location for AVIRIS or AVIRIS+MASTER that contains all flightlines
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
from shapely.geometry import shape, MultiPolygon
from shapely.ops import transform
from mpl_toolkits.basemap import Basemap
import mpl_toolkits.basemap as pyproj
from functools import partial

# Open shapefile
polygons = MultiPolygon([shape(pol['geometry']) for pol in fiona.open(polyLocation + '.shp')])
# project = partial(
#     pyproj.transform,
#     pyproj.Proj(polygonsOrg.crs),  # source coordinate system
#     pyproj.Proj(init='epsg:32611'))  # destination coordinate system
# polygons = shapely.ops.transform(project, polygonsOrg)  # apply projection
print len(polygons), "Polygons Found"

# Plot reference polygons
cm = plt.get_cmap('RdBu')
num_colours = len(polygons)
fig = plt.figure()
ax = fig.add_subplot(111)
minx, miny, maxx, maxy = polygons.bounds
w, h = maxx - minx, maxy - miny
ax.set_xlim(minx - 0.2 * w, maxx + 0.2 * w)
ax.set_ylim(miny - 0.2 * h, maxy + 0.2 * h)
ax.set_aspect(1)
patches = []
for idx, p in enumerate(polygons):
    colour = cm(1. * idx / num_colours)
    patches.append(PolygonPatch(p, fc=colour, ec='#555555', alpha=1., zorder=1))
ax.add_collection(PatchCollection(patches, match_original=True))
ax.set_xticks([])
ax.set_yticks([])
plt.title("Plant Species Reference Polygons")
plt.savefig('F:\\Image-To-Image-Registration\\AVIRIS\\test.png', alpha=True, dpi=300)
plt.show()
plt.clf()

# Find image files that spectra will be extracted from
flList = ['FL02']  # , 'FL03', 'FL04', 'FL05', 'FL06', 'FL07', 'FL08', 'FL09', 'FL10', 'FL11'
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
    imgFile = rasterio.open(i, 'r')
    imgProj = imgFile.affine  # Get projection
    print imgFile.read(1).shape

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
    rasterio.plot.show((imgFile, 137), cmap='magma')
    patches = []
    for idx, p in enumerate(polygons):
        colour = cm(1. * idx / num_colours)
        patches.append(PolygonPatch(p, fc=colour, ec='#555555', alpha=1., zorder=1))
    ax.add_collection(PatchCollection(patches, match_original=True))
    ax.set_xticks([])
    ax.set_yticks([])
    plt.title("Plant Species Reference Polygons")
    plt.savefig('F:\\Image-To-Image-Registration\\AVIRIS\\test.png', alpha=True, dpi=300)
    plt.show()

    polygonMask = rasterio.features.geometry_mask(polygons, out_shape=imgFile.shape, transform=imgFile.transform, all_touched=True)

    # Extract spectra from image
    # Note: cannot read entire AVIRIS file as it is too large, need to loop by band
    for bandNum in range(1, 2):  # bands are indexed from 1
        imgData = imgFile.read(bandNum)

        # create a masked numpy array
        # masked_data = np.ma.array(data=imgData, mask=polygonMask.astype(bool))
        # np.savetxt('F:\\Image-To-Image-Registration\\AVIRIS\\test.csv', masked_data, fmt="%s", delimiter=',')