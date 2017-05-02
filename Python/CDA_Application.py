# Created By: Susan Meerdink
# 4/24/17
# This python code apply canonical discriminant analysis eigenvectors to images
# --------------------------------------------------------------------------------------------------------------------

dirLocation = 'F:\\Image-To-Image-Registration\\AVIRIS\\'  # File location for AVIRIS or AVIRIS+MASTER that contains all flightlines
libLocation = 'F:\\Classification-Products\\3 - CDA Development\\'
dateTag = '140416'

import glob
import rasterio
import numpy as np

# Import CDA coefficients
libSpecCalFile = libLocation + dateTag + '_CDA_spectral_library_spectra.csv'
cda = np.loadtxt(libSpecCalFile, dtype=object, delimiter=',')  # Load in spectra - skips first line
cda = cda.astype(np.double)  # convert from string array to double array

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