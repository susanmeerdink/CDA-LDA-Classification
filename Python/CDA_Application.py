# Created By: Susan Meerdink
# 4/24/17
# This python code apply canonical discriminant analysis eigenvectors to images
# --------------------------------------------------------------------------------------------------------------------

dirLocation = 'F:\\Image-To-Image-Registration\\AVIRIS\\'  # File location for AVIRIS or AVIRIS+MASTER that contains all flightlines
libLocation = 'F:\\Classification-Products\\3 - CDA Development\\Combined Single Date\\'
dateTag = '140416'

import glob
import rasterio
import numpy as np
import matplotlib.pyplot as plt
import rasterio.plot
import gc

# Import CDA coefficients
libSpecCalFile = libLocation + dateTag + '_CDA_coefficients.csv'
cda = np.loadtxt(libSpecCalFile, dtype=object, delimiter=',')  # Load in spectra - skips first line
cda = cda.astype(np.double)  # convert from string array to double array

# Other Variables for analysis
bbl = np.array([0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1,
                1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0])

# Find image files
flList = ['FL02']  # 'FL02', 'FL03', 'FL04', 'FL05', 'FL06', 'FL07', 'FL08', 'FL09', 'FL10', 'FL11'
for fl in flList:  # loop through flightline files to find specific date
    imageLocation = dirLocation + fl + '\\6 - Spectral Correction Files\\*' + dateTag + '*'
    for name in glob.glob(imageLocation):  # Ignore header files
        if '.hdr' not in name:
            print 'Transforming to CDA Image from', name
            imgFile = rasterio.open(name, 'r')  # Open raster image
            shortName = name.split('\\')[-1]  # Get file name
            newImageLocation = dirLocation + fl + '\\7 - CDA Files\\' + shortName + '_CDA'
            outImgFile = rasterio.open(newImageLocation, 'w', width=imgFile.width, height=imgFile.height, count=22,
                                       transform=imgFile.transform, crs=imgFile.crs, driver='ENVI', dtype='float64')
            newData = np.empty([22, imgFile.height, imgFile.width])

            # Loop through columns of image and apply CDA coefficients
            # for col in range(0, imgFile.width):
            #     window = ((0, imgFile.height), (col, col + 1))
            #     origData = imgFile.read(window=window, masked=False, boundless=True)  # Returns spectra as bands x width x height
            #     origData = np.squeeze(origData)  # returns bands x height
            #     origData = np.delete(origData, np.where(bbl == 0)[0], 0)  # remove bad bands
            #     with np.errstate(all='ignore'):  # zero values in spectra will produce error, ignore them
            #         transData = np.log(origData)  # Transform image data using log (for normal distribution)
            #         transData[np.isneginf(transData)] = 0
            #     newData[:, :, col] = np.dot(cda, transData)

            # Loop through row of image and apply CDA coefficients
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
                origData = np.squeeze(origData)  # returns bands x width
                origData = np.delete(origData, np.where(bbl == 0)[0], 0)  # remove bad bands
                # with np.errstate(all='ignore'):  # zero values in spectra will produce error, ignore them
                #     transData = np.log(origData)  # Transform image data using log (for normal distribution)
                #     transData[np.isneginf(transData)] = 0
                # newData[:, row, :] = np.dot(cda.T, transData)
                newData[:, row, :] = np.dot(cda.T, origData)

            # fig, axMap = plt.subplots(num=None, figsize=(4, 3), dpi=300, facecolor='w', edgecolor='k')
            # plt.imshow(newData[1,:,:])
            # plt.show()

            # Save new data into raster file
            print 'Writing data for ' + newImageLocation
            outImgFile.write(newData)  # write and save new CDA image

            print 'Completed processing for ' + newImageLocation
            outImgFile.close()
