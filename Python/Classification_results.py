# Created By: Susan Meerdink
# 5/9/17
# This python code reads in LDA classifications results and calculates results
# --------------------------------------------------------------------------------------------------------------------

# Variable Names
dirLocation = 'F:\\Image-To-Image-Registration\\AVIRIS\\'  # File location for AVIRIS that contains all flightlines
outLocation = 'F:\\Classification-Products\\4 - LDA Classification\\Combined Single Date\\'
polyLocation = 'F:\\Classification-Products\\2017_05_01_reference_polygons_withMeta'
metaLocation = 'F:\\Classification-Products\\2017_05_01_reference_polygons_catalog.csv'
dateTag = '140416'

# Import Modules and Libraries
import rasterio.plot
from rasterio import features
import glob
import numpy as np
import fiona
from shapely.geometry import shape, MultiPolygon
import pyproj
from functools import partial
from shapely.ops import transform
import gc

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

# Misc Variables
classNames = ('ADFA','AGRES','ARCA-SALE','ARGL','BAPI','BRNI','CECU','CEME','CESP','CISP','ERFA','EUSP','IRGR','MAGF',
              'PEAM','PISA','PLRA','QUAG','QUDO','ROCK','SOIL','UMCA','URBAN')

# Variables to hold the entire spectral library, calibration library, and validation library with metadata
spectralLibData = np.empty([0, 2])
spectralLibName = np.empty([0, 5])
spectralLibMeta = np.empty([0, len(headers) + 6])

# Find image files
flList = ['FL02']  # 'FL02', 'FL03', 'FL04', 'FL05', 'FL06', 'FL07', 'FL08', 'FL09', 'FL10', 'FL11'
for fl in flList:  # loop through flightline files to find specific date
    imageLocation = dirLocation + fl + '\\8 - LDA Classification Files\\*' + dateTag + '*'
    for name in glob.glob(imageLocation):  # Ignore header files
        if '.hdr' not in name and '.xml' not in name:
            print 'Calculation of LDA classification results for ', name
            imgFile = rasterio.open(name, 'r')  # Open raster image
            shortName = name.split('\\')[-1]  # Get file name

            for idx in range(0, len(polygons)):  # Loop through polygons
                gc.collect()
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
                            inName = [fl, dateTag, polyName, x, y]
                            inMeta = np.hstack((inName, metadata[np.where(polyIndex == polyName)[0][0], :], pixel))
                            inData = np.hstack((int(metadata[np.where(polyIndex == polyName)[0][0], 10]), int(pixel)))
                            spectralLibData = np.vstack((spectralLibData, inData))
                            spectralLibName = np.vstack((spectralLibName, inName))
                            spectralLibMeta = np.vstack((spectralLibMeta, inMeta))

            # Output Classification Results
            headerOut = 'Flightline, Date, PolygonName, X, Y,' + ','.join(headers) + ', Classification'
            fileOut = file((outLocation + dateTag + '_LDA_Classification_Results_by_Pixel.csv'), 'wb')
            np.savetxt(fileOut, spectralLibMeta, header=headerOut, fmt='%s', delimiter=',')

            # Create Classification Matrix
            classMatrix = np.zeros([23, 23])
            for dom in range(1, 24):
                oneDom = spectralLibData[np.where(spectralLibData[:, 0] == dom)[0], :]
                if oneDom.shape[0] > 0:
                    for i in range(0, oneDom.shape[0]):
                        classMatrix[dom-1, int(oneDom[i,1])-1] = classMatrix[dom-1, int(oneDom[i,1])-1] + 1

            # Summing up Values#
            rowTotals = np.sum(classMatrix, axis=1)  # sum up each row
            columnTotals = np.sum(classMatrix, axis=0)  # sum up each column
            pixelTotal = np.sum(rowTotals)  # Get total number of pixels

            # Variables for Classification User's and Producer's Accuracy#
            count = 0  # Set up value to loop through
            producer = np.zeros(23)  # Create variable that will hold producer's accuracy
            user = np.zeros(23)  # Create variable that will hold producer's accuracy
            totalCorrect = 0  # Variable that will sum up total number of correctly classified pixels

            # Classification User's and Producer's Accuracy
            while count < classMatrix.shape[0]:  # While we haven't reached end of classification array
                classCorrect = classMatrix[count, count]  # get the total number of pixels for class that have been classified correctly
                if classCorrect == 0:
                    producer[count] = 0
                    user[count] = 0
                else:
                    producer[count] = classCorrect / columnTotals[count]  # Divide the correct number of pixels by the total number of actual pixels present
                    user[count] = classCorrect / rowTotals[count]  # Divide the correct number of pixels by the total number of actual pixels present

                totalCorrect = totalCorrect + classCorrect  # add current correctly classified pixels to total correct pixels
                count = count + 1  # increase the counter

            # Overall Accuracy
            overall = totalCorrect / pixelTotal  # divide the total number of correctly classified pixels by the total number of pixels

            # Kappa
            topKappa = (totalCorrect * pixelTotal) - np.sum(rowTotals, axis=0)  # Multiply correctly classified pixels with total number of pixels subtract row totals
            bottomKappa = (pixelTotal * pixelTotal) - np.sum(rowTotals, axis=0)  # Multiply total number of pixels with total number of pixels subtract row totals
            kappa = topKappa / bottomKappa  # Get Kappa statistic by dividing top and bottom

            # Write Results
            outputFile = open((outLocation + dateTag + '_LDA_Classification_Matrix.csv'), 'w')  # open the output file for writing
            outputFile.write(',' + ','.join(classNames) + ',Row Total\n')  # Write first line with header
            loopRange = range(0, classMatrix.shape[0])  # Set range of values to loop through
            for value in loopRange:
                outputFile.write(classNames[value] + ',')  # Write classification name
                outputFile.write(','.join(map(str,classMatrix[value, :])))  # Write out classification array
                outputFile.write(',' + str(rowTotals[value]) + '\n')  # Finish row writing with row total

            outputFile.write('Column Total,' + ','.join(map(str, columnTotals)) + ',' + str(pixelTotal) + '\n')  # Write out column totals
            outputFile.write('\n')  # Add another blank row for visual
            outputFile.write('Producers Accuracy,' + ','.join(map(str, producer)) + '\n')  # Write out producer's accuracy
            outputFile.write('Users Accuracy,' + ','.join(map(str, user)) + '\n')  # Write out user's accuracy
            outputFile.write('Number of Correctly Classified Pixels,' + str(totalCorrect) + '\n')  # Write total number of correct pixels
            outputFile.write('Overall Accuracy,' + str(overall) + '\n')  # Write overall accuracy
            outputFile.write('Kappa,' + str(kappa) + '\n')  # Write out kappa statistic

            print 'Finished Classification Accuracy Results for: ', name
            outputFile.close()

