# Susan Meerdink
# 3/18/17
# This file reads in a envi spectral library and splits it into independent validation
# and model development. Then splits the model development into training/validation
############INPUTS######################################
specLibFile = 'F:\\Classification-Products\\FL03\\1 - Spectral Library\\f140829_AVIRIS_spectral_library' #The envi spectral library
specMetaFile = 'F:\\Classification-Products\\FL03\\1 - Spectral Library\\f140829_AVIRIS_spectral_library.csv' #The envi spectral library

############ENDINPUTS####################################

import numpy as np
import spectral.io.envi as envi
import spectral
import random
import re
from sets import Set
#from sklearn.discriminant_analysis import LinearDiscriminantAnalysis

specLib = envi.open(specLibFile+'.hdr') #Open envi spectral library
spectra = specLib.spectra #pull out spectra
metadata = np.loadtxt(specMetaFile,dtype=object,delimiter = ',')#Load in metadata
headers = metadata[0,:] #save headers separate of metadata
headers = np.char.strip(headers.astype(str)) #remove whitespace from headers
#sortIndex = np.where(headers == 'Dominant')[0][0] #pull out what class to stratify sampling by
polyIndex = np.where(headers == 'Polygon_ID')[0][0] #pull out what class to stratify sampling by
metadata = np.delete(metadata,0,0) #remove the headers
#species = np.unique(metadata[:,sortIndex])#Get a list of classes/species in library
polygons = np.unique(metadata[:,polyIndex])#Get list of polygons in library

#Create new output files
calLibSpec = file(specLibFile + '_spectra_calibration.csv','a')
calLibMeta = file(specLibFile + '_metadata_calibration.csv','a')
valLibSpec = file(specLibFile + '_spectra_validation.csv','a')
valLibMeta = file(specLibFile + '_metadata_validation.csv','a')
calLibMeta.write(','.join(headers) + '\n')
valLibMeta.write(','.join(headers) + '\n')

##Loop through polygons
for p in polygons:
    print('Processing: ' + p)
    spPixels = np.where(metadata[:,polyIndex] == p) # find pixels for this polygon
    spMeta = metadata[spPixels[0],:] #pull out pixels metadata
    spSpec = spectra[spPixels[0],:] #pull out pixels spectra
    print('Number of Pixels: ' + str(len(spPixels[0])))
    valIndex = []

    #Using Proportional Limit of 50% for smaller polygons and for
    #Absolute Limit of 10 spectra for larger polygons (Roth et al. 2012)
    propLimit = 0.5
    absoLimit = 10

    #Separate into validation/calibration
    if propLimit* len(spPixels[0]) < 10: #If it is a small polygon, use proportional limit
        valIndex = random.sample(xrange(0, len(spPixels[0])-1),int(round(propLimit* len(spPixels[0]))))
    else: #If it is a large polygon use absolute limit
        valIndex = random.sample(xrange(0, len(spPixels[0])-1),absoLimit)
    fullIndex = range(0,len(spMeta)-1)
    calIndex = list(Set(fullIndex).difference(valIndex)) 
    pValMeta = spMeta[valIndex,:]
    pValSpec = np.column_stack((pValMeta[:,0],spSpec[valIndex,:])).T
    pCalMeta = spMeta[calIndex,:]
    pCalSpec = np.column_stack((pCalMeta[:,0],spSpec[calIndex,:])).T

    #store into new files
    np.savetxt(calLibSpec,pCalSpec,fmt = "%s", delimiter = ',')
    np.savetxt(calLibMeta,pCalMeta,fmt = "%s", delimiter = ',')
    np.savetxt(valLibSpec,pValSpec,fmt = "%s", delimiter = ',')
    np.savetxt(valLibMeta,pValMeta,fmt = "%s", delimiter = ',')

calLibSpec.close()
calLibMeta.close()
valLibSpec.close()
valLibMeta.close()
        
#Loop through classes/species
##for sp in species:
##    print('Processing: ' + sp)
##    spPixels = np.where(metadata[:,sortIndex] == sp) #find pixels for this species
##    spMeta = metadata[spPixels[0],:] #pull out pixels metadata
##    spSpec = spectra[spPixels[0],:] #pull out pixels spectra
##    spPolyList = np.unique(spMeta[:,polyIndex]) #pull out unique polygons for this species
##    print('Number of Polygons: ' + str(len(spPolyList)))
##    valIndex = []
##
##    #Split up into validation and calibration
##    if len(spPolyList) > 2: #if there is more than two polygon
##        numValPoly = int(round(0.2 * len(spPolyList))) #calculate the number of polygons for independent validation
##        randIndex = random.sample(xrange(1, len(spPolyList)),numValPoly)
##        for i in range(0,len(randIndex)):
##            valIndex.extend(np.where(spMeta[:,polyIndex] == spPolyList[randIndex[i]])[0])
##        fullIndex = range(0,len(spMeta))
##        calIndex = list(Set(fullIndex).difference(valIndex))
##        spValMeta = spMeta[valIndex,:]
##        spValSpec = spSpec[valIndex,:]
##        spCalMeta = spMeta[calIndex,:]
##        spCalSpec = spSpec[calIndex,:]
##
##    elif len(spPolyList) == 2: #if there is two polygons
##        numValPoly = 1
##        randIndex = random.sample(xrange(1, 2), 1)
##        pullIndex = np.where(spMeta[:,polyIndex] == spPolyList[randIndex])[0]
##        fullIndex = range(0,len(spMeta))
##        calIndex = list(Set(fullIndex).difference(valIndex))
##        spValMeta = spMeta[valIndex,:]
##        spValSpec = spSpec[valIndex,:]
##        spCalMeta = spMeta[calIndex,:]
##        spCalSpec = spSpec[calIndex,:]
##                                  
##    else: #if there is one polygon
##        numValPoly = 0
##        spValMeta = 0
##        spValSpec = 0
##        spCalMeta = spMeta[range(0,len(spMeta)),:]
##        spCalSpec = spSpec[range(0,len(spMeta)),:]

    

    

