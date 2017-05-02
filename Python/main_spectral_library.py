# Created By: Susan Meerdink
# 3/27/17
# This python code is the main script for development of spectral libraries
# --------------------------------------------------------------------------------------------------------------------
# Inputs
polyLocation = 'F:\\Classification-Products\\2017_05_01_reference_polygons_withMeta'
metaLocation = 'F:\\Classification-Products\\2017_05_01_reference_polygons_catalog.csv'
dirLocation = 'F:\\Image-To-Image-Registration\\AVIRIS\\*'  # File location for AVIRIS or AVIRIS+MASTER that contains all flightlines
outLoc = 'F:\\Classification-Products\\1 - Spectral Library\\'  # Ouput file location for figures generated from analysis
dateTag = '140416'

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
import pyproj
from functools import partial
import random
from sets import Set
from statsmodels import robust

# Open shapefile and transform to appropriate coordinate system
print "Opening reference polygons from " + polyLocation
polyCRS = fiona.open(polyLocation + '.shp', "r")  # Open polygon to get crs information for transformation
polyOriginal = MultiPolygon([shape(pol['geometry']) for pol in fiona.open(polyLocation + '.shp')])  # Open polygon as a MultiPolygon for processing purposes
project = partial(  # Define function for projection process
    pyproj.transform,
    pyproj.Proj(polyCRS.crs),  # source coordinate system
    pyproj.Proj(init='epsg:32611'))  # destination coordinate system, UTM Zone 11 WGS 84
polygons = transform(project, polyOriginal)  # apply projection
print len(polygons), "Reference Polygons Found"

# Load Metadata
metadata = np.loadtxt(metaLocation, dtype=object, delimiter=',')  # Load in metadata
headers = metadata[0, :]  # save headers separate of metadata
headers = np.char.strip(headers.astype(str))  # remove whitespace from headers
metadata = np.delete(metadata, 0, 0)  # remove the headers
polyIndex = metadata[:, (np.where(headers == 'PolyID')[0][0])]  # pull out what class to stratify sampling by

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

# Variables to hold the entire spectral library, calibration library, and validation library with metadata
spectralLibData = np.empty([0, 224])
spectralLibName = np.empty([0, 5])
spectralLibMeta = np.empty([0, len(headers) + 5])
valMeta = np.empty([0, len(headers) + 5])
calMeta = np.empty([0, len(headers) + 5])
valSpec = np.empty([0, 229])
calSpec = np.empty([0, 229])

# Find image files and extract spectra using reference polygons
flList = ['FL02', 'FL03', 'FL04', 'FL05', 'FL06', 'FL07', 'FL08', 'FL09', 'FL10', 'FL11']  # 'FL02', 'FL03', 'FL04', 'FL05', 'FL06', 'FL07', 'FL08', 'FL09', 'FL10', 'FL11'
offset = 0  # counter used to pull out  validation and calibration
for fl in flList:  # loop through flightline files to find specific date
    imageLocation = dirLocation + fl + '\\6 - Spectral Correction Files\\*' + dateTag + '*'
    for name in glob.glob(imageLocation):  # Ignore header files
        if '.hdr' not in name:
            print 'Extracting spectra from', name
            imgFile = rasterio.open(name, 'r')  # Open raster image
            shortName = name.split('\\')[-1]  # Get file name
            spectraCount = 0

            # Plot reference polygons on top of current image
            fig, axMap = plt.subplots(num=None, figsize=(4, 3), dpi=300, facecolor='w', edgecolor='k')
            rasterio.plot.show((imgFile, 137), ax=axMap, cmap='hot', clim=(1, 10000))  # Show band 137 of current image
            patches = []
            for idx, p in enumerate(polygons):
                patches.append(PolygonPatch(p, fc='#FFFFFF', ec='#FFFFFF', alpha=1., zorder=1))
            axMap.add_collection(PatchCollection(patches, match_original=True))
            plt.xticks(fontsize=6)
            plt.yticks(fontsize=6)
            plt.title("Plant Species Reference Polygons with \n" + shortName, fontsize=6)
            plt.savefig(outLoc + 'Figures\\' + shortName + '_reference_poly.png', alpha=True, dpi=300)
            plt.clf()

            for idx in range(0, len(polygons)):  # Loop through polygons
                polyIn = polygons[idx]
                polyName = polyCRS[idx]['properties']['PolyID']
                pixelCount = 0

                # Create Mask that has 1 for locations with polygons and 0 for non polygon locations
                polygonMask = rasterio.features.rasterize([(polyIn, 1)], out_shape=imgFile.shape,
                                                          transform=imgFile.transform, all_touched=False)
                test = np.count_nonzero(polygonMask)  # Get the number of elements that are not zero
                if test > 0:  # If there is data for this polygon assign the data
                    indices = np.nonzero(polygonMask)
                    for i in range(0, len(indices[0])):  #
                        x = indices[0][i]
                        y = indices[1][i]
                        window = ((x, x + 1), (y, y + 1))
                        data = imgFile.read(window=window, masked=False, boundless=True)  # Extract spectra from image
                        pixel = np.transpose(data[:, 0, 0])
                        if any(pixel):  # If there are non zero values save them to spectral library
                            pixelCount += 1  # How many pixels are in this polygon
                            spectraCount += 1  # How many spectra were collected from flightline
                            inName = [fl, dateTag, polyName, x, y]
                            inMeta = np.hstack((inName, metadata[np.where(polyIndex == polyName)[0][0], :]))
                            spectralLibData = np.vstack((spectralLibData, pixel))
                            spectralLibName = np.vstack((spectralLibName, inName))
                            spectralLibMeta = np.vstack((spectralLibMeta, inMeta))

                    if pixelCount > 0:  # Split spectra into training and validation libraries
                        # Pull out bad spectra
                        # print spectralLibMeta[offset+1, 2]
                        # fullIndex = range(offset, pixelCount + offset)
                        # mad = robust.mad(spectralLibData[fullIndex, :], c=1)  # http://www.statsmodels.org/stable/generated/statsmodels.robust.scale.mad.html#statsmodels.robust.scale.mad
                        # diff = np.zeros(spectralLibData[fullIndex,:].shape)
                        # index = np.where(spectralLibData[fullIndex, :] > (mad * 3) + np.median(spectralLibData[fullIndex, :], axis=0))
                        # diff[index[0], index[1]] = 1
                        # index = np.where(spectralLibData[fullIndex, :] < np.median(spectralLibData[fullIndex, :], axis=0) - (mad * 3))
                        # diff[index[0], index[1]] = 1
                        # sumArray = np.sum(diff, axis=1)
                        # outlierIndex = np.where(sumArray > 112)
                        # outlierIndex = [x + offset for x in outlierIndex]

                        # spectralLibData = np.delete(spectralLibData, outlierIndex[0], axis=0)
                        # spectralLibMeta = np.delete(spectralLibMeta, outlierIndex[0], axis=0)
                        # spectralLibName = np.delete(spectralLibName, outlierIndex[0], axis=0)
                        # pixelCount = pixelCount - len(outlierIndex[0])
                        # spectraCount = spectraCount - len(outlierIndex[0])

                        # Using Proportional Limit of 50% for smaller polygons and for
                        # Absolute Limit of 10 spectra for larger polygons (Roth et al. 2012)
                        propLimit = 0.5
                        absoLimit = 10

                        # Separate into validation/calibration
                        if propLimit * pixelCount < 10:  # If it is a small polygon, use proportional limit
                            calIndex = random.sample(xrange(0, pixelCount), int(round(propLimit * pixelCount)))
                        else:  # If it is a large polygon use absolute limit
                            calIndex = random.sample(xrange(0, pixelCount), absoLimit)
                        fullIndex = range(0, pixelCount)
                        valIndex = list(Set(fullIndex).difference(calIndex))
                        valIndex = [x + offset for x in valIndex]
                        calIndex = [x + offset for x in calIndex]
                        valMeta = np.vstack((valMeta, spectralLibMeta[valIndex, :]))
                        calMeta = np.vstack((calMeta, spectralLibMeta[calIndex, :]))
                        inValSpec = np.hstack((spectralLibName[valIndex, :], spectralLibData[valIndex, :]))
                        inCalSpec = np.hstack((spectralLibName[calIndex, :], spectralLibData[calIndex, :]))
                        valSpec = np.vstack((valSpec, inValSpec))
                        calSpec = np.vstack((calSpec, inCalSpec))
                    offset = offset + pixelCount
                else:  # If there is no data for this polygon skip to the next one
                    continue
            print 'Done extracting ', spectraCount, ' spectra'

# Create new output files
fileOutSpec = file((outLoc + 'Combined Single Date\\' + dateTag + '_spectral_library_spectra.csv'), 'wb')
fileOutMeta = file((outLoc + 'Combined Single Date\\' + dateTag + '_spectral_library_metadata.csv'), 'wb')
calLibSpec = file(outLoc + 'Combined Single Date\\' + dateTag + '_spectral_library_calibration_spectra.csv', 'wb')
calLibMeta = file(outLoc + 'Combined Single Date\\' + dateTag + '_spectral_library_calibration_metadata.csv', 'wb')
valLibSpec = file(outLoc + 'Combined Single Date\\' + dateTag + '_spectral_library_validation_spectra.csv', 'wb')
valLibMeta = file(outLoc + 'Combined Single Date\\' + dateTag + '_spectral_library_validation_metadata.csv', 'wb')

# Get together data and headers
headerOutSpec = 'Flightline, Date, PolygonName, X, Y,' + ','.join(map(str, imgFile.indexes))
headerOutMeta = 'Flightline, Date, PolygonName, X, Y,' + ','.join(headers)
allSpec = np.hstack((spectralLibName, spectralLibData))
allMeta = spectralLibMeta

# Save spectral libraries from all flightlines
np.savetxt(fileOutSpec, allSpec, header=headerOutSpec, fmt='%s', delimiter=",")
np.savetxt(fileOutMeta, allMeta, header=headerOutMeta, fmt='%s', delimiter=",")
np.savetxt(calLibSpec, calSpec, header=headerOutSpec, fmt="%s", delimiter=',')
np.savetxt(calLibMeta, calMeta, header=headerOutMeta, fmt="%s", delimiter=',')
np.savetxt(valLibSpec, valSpec, header=headerOutSpec, fmt="%s", delimiter=',')
np.savetxt(valLibMeta, valMeta, header=headerOutMeta, fmt="%s", delimiter=',')

# Close files
fileOutSpec.close()
fileOutMeta.close()
calLibSpec.close()
calLibMeta.close()
valLibSpec.close()
valLibMeta.close()
