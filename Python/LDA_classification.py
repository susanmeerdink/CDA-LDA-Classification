# Created By: Susan Meerdink
# 5/9/17
# This python code classifies CDA images into plant species
# --------------------------------------------------------------------------------------------------------------------

#Variable Names
dirLocation = 'F:\\Image-To-Image-Registration\\AVIRIS\\'  # File location for AVIRIS that contains all flightlines
libLocation = 'F:\\Classification-Products\\3 - CDA Development\\Combined Single Date\\'  # File location for CDA library
outLocation = 'F:\\Classification-Products\\4 - LDA Classification\\Combined Single Date\\'
dateTag = '140416'

# Import Modules/Libraries
import glob
import rasterio
import numpy as np
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn import linear_model
import matplotlib.pyplot as plt
import gc

# Read in metadata files
libMetaCalFile = libLocation + dateTag + '_CDA_spectral_library_calibration_metadata.csv'
libMetaValFile = libLocation + dateTag + '_CDA_spectral_library_validation_metadata.csv'
metaCal = np.loadtxt(libMetaCalFile, dtype=object, delimiter=',')
metaVal = np.loadtxt(libMetaValFile, dtype=object, delimiter=',')

# Read in calibration spectral files
libSpecCalFile = libLocation + dateTag + '_CDA_spectral_library_calibration_spectra.csv'
spectraCal = np.loadtxt(libSpecCalFile, dtype=object, delimiter=',')  # Load in spectra - skips first line
metaSpecCal = spectraCal[:, 0:5]  # save first 5 columns of spectra separately
spectraCal = np.delete(spectraCal, [0, 1, 2, 3, 4], 1)  # remove the 5 columns of metadata in spectra
spectraCal = spectraCal.astype(np.double)  # convert from string array to double array

# Read in validation spectral files
libSpecValFile = libLocation + dateTag + '_CDA_spectral_library_validation_spectra.csv'
spectraVal = np.loadtxt(libSpecValFile, dtype=object, delimiter=',')  # Load in spectra - skips first line
metaSpecVal = spectraVal[:, 0:5]  # save first 5 columns of spectra separately
spectraVal = np.delete(spectraVal, [0, 1, 2, 3, 4], 1)  # remove the 5 columns of metadata in spectra
spectraVal = spectraVal.astype(np.double)  # convert from string array to double array

# LDA classification training
clf = LinearDiscriminantAnalysis()  # http://scikit-learn.org/stable/modules/generated/sklearn.discriminant_analysis.LinearDiscriminantAnalysis.html#sklearn.discriminant_analysis.LinearDiscriminantAnalysis
clf.fit(spectraCal, metaCal[:, 15].astype(np.int))
predict = clf.predict(spectraVal)
regr = linear_model.LinearRegression()
x = metaVal[:, 15].astype(np.int).reshape(len(metaVal[:, 15]), 1)
y = predict.reshape(len(predict), 1)
linearResults = regr.fit(x, y)  # http://scikit-learn.org/stable/modules/generated/sklearn.linear_model.LinearRegression.html
plt.scatter(metaVal[:, 15].astype(np.int), predict)
plt.plot(metaVal[:, 15].astype(np.int), regr.predict(x))
plt.ylabel('Predicted')
plt.xlabel('Observed')
plt.title('LDA ' + dateTag)
plt.ylim([0, 25])
plt.xlim([0, 25])
plt.savefig(outLocation + dateTag + '_LDA_Pred_v_Obs.png')
ldaResults = np.empty([1, 4])  # Array for intercept, slope, r2, and rmse
ldaResults[0, 0] = linearResults.intercept_  # Intercept
ldaResults[0, 1] = linearResults.coef_  # Slope
ldaResults[0, 2] = linearResults.score(x, predict.reshape(len(predict), 1))  # R2
ldaResults[0, 3] = np.mean(regr.predict(x - predict.reshape(len(predict), 1)) ** 2)  # RMSE

# Find image files
flList = ['FL02']  # 'FL02', 'FL03', 'FL04', 'FL05', 'FL06', 'FL07', 'FL08', 'FL09', 'FL10', 'FL11'
for fl in flList:  # loop through flightline files to find specific date
    imageLocation = dirLocation + fl + '\\7 - CDA Files\\*' + dateTag + '*'
    for name in glob.glob(imageLocation):  # Ignore header files
        if '.hdr' not in name and '.xml' not in name:
            print 'LDA Classification of CDA Image ', name
            imgFile = rasterio.open(name, 'r')  # Open raster image
            shortName = name.split('\\')[-1]  # Get file name
            newImageLocation = dirLocation + fl + '\\8 - LDA Classification Files\\' + shortName + '_LDA'
            outImgFile = rasterio.open(newImageLocation, 'w', width=imgFile.width, height=imgFile.height, count=1,
                                       transform=imgFile.transform, crs=imgFile.crs, driver='ENVI', dtype='float64')
            newData = np.empty([imgFile.height, imgFile.width])

            # Loop through row of image
            for row in range(0, imgFile.height):
                gc.collect()
                if row == int(imgFile.height*0.25):
                    print '25% Completed'
                if row == int(imgFile.height*0.50):
                    print '50% Completed'
                if row == int(imgFile.height*0.75):
                    print '75% Completed'
                window = ((row, row + 1), (0, imgFile.width))
                origData = imgFile.read(window=window, masked=False, boundless=True)  # Returns spectra as bands x width x height
                origData = np.squeeze(origData).T
                if np.any(np.isnan(origData)) == 1:
                    origData[np.argwhere(np.isnan(origData))] = 0
                unclassIndex = np.all(origData == 0, 1)
                classData = clf.predict(origData)
                classData[np.where(unclassIndex == True)] = 0
                newData[row, :] = classData

            # Save new data into raster file
            print 'Writing data for ' + newImageLocation
            outImgFile.write(newData, 1)  # write and save new CDA image
            print 'Completed processing for ' + newImageLocation
            outImgFile.close()

# Save linear regression results from Canonical Variables to Text File
outCDA = file(outLocation + dateTag + '_LDA_classification_results.csv', 'wb')
headerOut = 'Intercept, Slope, R^2, RMSE'
np.savetxt(outCDA, ldaResults, header=headerOut, fmt='%s', delimiter=',')
outCDA.close()