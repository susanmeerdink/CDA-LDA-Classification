# Susan Meerdink
# updated 12/8/2015
# last updated 3/18/2016
# ############INPUTS######################################
#MUST add semi colon (;) to end of file!
OrigROI = 'R:\\users\\susan.meerdink\\Dropbox\\AAG_2016_Research\\ROIs\\20130411_AVIRIS_diff.txt'
#OrigROI = 'R:\\users\\susan.meerdink\\Dropbox\\AAG_2016_Research\\ROIs\\20140829_AVIRIS_same.txt'
#OrigROI = 'R:\\users\\susan.meerdink\\Dropbox\\AAG_2016_Research\\ROIs\\20130411_AVIRIS&MASTER_diff.txt'
#OrigROI = 'R:\\users\\susan.meerdink\\Dropbox\\AAG_2016_Research\\ROIs\\20130411_AVIRIS&MASTER_same.txt'

OutputROI = 'R:\\users\\susan.meerdink\\Dropbox\\AAG_2016_Research\\ROIs\\20130411_AVIRIS_diff_OUTPUT.csv'
#OutputROI = 'R:\\users\\susan.meerdink\\Dropbox\\AAG_2016_Research\\ROIs\\20140829_AVIRIS_same_OUTPUT.csv'
#OutputROI = 'R:\\users\\susan.meerdink\\Dropbox\\AAG_2016_Research\\ROIs\\20140829_AVIRIS&MASTER_diff_OUTPUT.csv'
#OutputROI = 'R:\\users\\susan.meerdink\\Dropbox\\AAG_2016_Research\\ROIs\\20140606_AVIRIS&MASTER_same_OUTPUT.csv'

MasterMetaData = 'R:\\users\\susan.meerdink\\Dropbox\\AAG_2016_Research\\Spectral Libraries\\SBFR_all432polys_2013_metadata.csv'

############ENDINPUTS####################################

import numpy as np
#Set Output file for ROIs
outputFile = open(OutputROI,'w')
inputFile = open(OrigROI,'r')

meta = [] #Variable to hold metadata
i = 0 #counter
#Loop through file that contains all the metadata fields for each polygon and store the metadata so that it can be assigned to the new library
for polygon in open(MasterMetaData):
    text = polygon.split(',')#Split the line into array based on comma
    string = text[-1]#Get last element of list
    top = string.find('\n') #Find the newline character
    text[-1] = string[0:top] #Remove the newline character from the last element
    meta.append(text) #Add this line (polygon info) to the meta data variable 
    i = i +1 #Advance counter

#Create headers for the output file
outputFile.write('Name,# of Pts,')
outputFile.write( ', '.join(meta[0]))
outputFile.write(',ID,X,Y, Map X, Map Y, Lat, Long, Classification Code\n')

polygonIDList = []
ptList = []
roiData = []
line = inputFile.readline()
while ';' in line:
    #Finding and storing polygon name
    bottomPoly = line.find('=')
    if bottomPoly > -1:
        topPoly = line.find(')')
        polygonID = line[bottomPoly+1:topPoly]
        polygonIDList.append(polygonID)
        #outputFile.write(polygonID + ',')

    #Finding and Storing number of points
    bottomPt = line.find('npts:')
    if bottomPt > -1:
        inLine = line[bottomPt+5:bottomPt+8]
        pt = inLine.strip()
        #outputFile.write(pt + '\n')
        ptList.append(pt)

    line = inputFile.readline()

k = 0
number = range(0,i)#create a list with numbers ranging from 0 to the length of the polygon metadata (first file that was looped through)
while ';' not in line:
    #If it's a blank line, then go on to the next polygon ID
    if line.isspace() == 1:
        k = k +1
        line = inputFile.readline()

    else:
        #Write out polygon name and points
        outputFile.write(polygonIDList[k] + ',' + ptList[k]+ ',')

        #Write out metadata fields
        for j in number: #Loop through polygon metadata
            if polygonIDList[k] == meta[j][1]: #Check to see if the pixel belongs to the metadata polygon, if it does
                outputFile.write(', '.join(meta[j])+ ',')#
                break
            j = j +1

        #Write out line
        #string = filter(None,line.split('\t'))
        string = filter(None,line.split(' '))
        outputFile.write(','.join(string))
        line = inputFile.readline()
    
outputFile.close()
inputFile.close()
    

#########################END#####################################################
    
