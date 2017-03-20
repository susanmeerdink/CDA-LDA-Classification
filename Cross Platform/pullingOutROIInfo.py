# Extracting ROI data
# Susan Meerdink
# Created 12/8/2015
# This code reads through a spectral library metdata file created with ENVI and formats it into .csv.
# Each row of the csv is one item in the spectral library. Does not include RGB value of ROI.
# Does include Name, # of Points, spectral information, and other information included in ENVI output
# (such as X, Y, MapX, MapY, Lat, Lon).
# This code takes two inputs:
# 1. The original spectral library from ENVI
# 2. Output file name

# ############INPUTS######################################
OrigROI = 'I:\\Classification-Products\FL03\\3 - Classification Results\\f140829_class_results.txt'

OutputROI = 'I:\\Classification-Products\FL03\\3 - Classification Results\\f140829_class_results_OUTPUT.csv'

############ENDINPUTS####################################

import numpy as np

outputFile = open(OutputROI,'w') #Create the output file for the formatted ROIs
inputFile = open(OrigROI,'r') #Open the original ROI file 

## Read through top portion of ENVI Metadata to pull out Name, Number of points, and header
polygonIDList = [] #Will hold the polygon list
ptList = [] #Will hold the number of points for each polygon
k = 0 #Use this to loop through polygon list

#line = inputFile.readline() #Sets the inital line to be read in
for line in inputFile:
    #print(line)
    if ';' in line: #Pull out header info
        #Finding and storing polygon name
        bottomPoly = line.find('name:')
        if bottomPoly > -1:
            topPoly = line.find('\n')
            polygonID = line[bottomPoly+6:topPoly]
            polygonIDList.append(polygonID)

        #Finding and Storing number of points
        bottomPt = line.find('npts:')
        if bottomPt > -1:
            inLine = line[bottomPt+5:bottomPt+9]
            pt = inLine.strip()
            ptList.append(pt)

        #Finding and Writing header for output file
        header = line.find(' ID ')
        if header > -1:
            headerLine = line.split('  ')
            headerLineIn = [x for x in headerLine if x]
            outputFile.write('Name,nPts,')
            outputFile.write(', '.join(headerLineIn[1:len(headerLineIn)]))#Write header line to file (leaving out semi colon as first item)

    else: #Pull out Spectral Information
        #If it's a blank line, then go on to the next polygon ID
        if line.isspace() == 1:
            k = k +1

        else:            
            #Write out polygon name and points
            outputFile.write(polygonIDList[k] + ',' + ptList[k]+ ',')

            #Write out line
            string = filter(None,line.split(' '))
            outputFile.write(','.join(string))
    
outputFile.close()
inputFile.close()
print('Formatting Complete')

#########################END#####################################################
    
