# Susan Meerdink
# last updated 12/1/2015
# This file reads in a spectral library metadata folder that is not filled out. This spectral library will have been generated using Viper tools 2.0 'build spectral library from rois'.
# It takes in a master file that contains all the information for each polygon and then assigns that polygon information to each pixel spectral in the unfinished spectral library.
############INPUTS######################################
UncompleteSpecLib = 'R:\\users\\susan.meerdink\\Dropbox\\AAG_2016_Research\\Spectral Libraries\\20140829_Spectral_Library_AVIRIS&MASTER_uncompleted.csv' #The library that is missing metadata
MasterMetaData = 'R:\\users\\susan.meerdink\\Dropbox\\AAG_2016_Research\\Spectral Libraries\\SBFR_all432polys_2013_metadata.csv' #This excel sheet contains the metadata
OutputSpecLib = 'R:\\users\\susan.meerdink\\Dropbox\\AAG_2016_Research\\Spectral Libraries\\20140829_Spectral_Library_AVIRIS&MASTER.csv' #This will be the final library WITH metadata
############ENDINPUTS####################################
############OPTIONAL INPUTS##############################
resolution = '18'
date = '8/29/2014'
sensor = 'AVIRIS & MASTER'
flightline = '2 & 3'
############END OF OPTIONAL INPUTS#######################
import numpy as np
#Set Output file for spectral library
outputFile = open(OutputSpecLib,'w')

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
outputFile.write('Name,')
outputFile.write( ', '.join(meta[0]))
outputFile.write(',Resolution_m,'+'Date,'+'Sensor,'+'FlightLine\n')


sliMeta = [] #Variable to hold metadata file
#Loop through spectral library that does NOT contain completed metadata fields
for pixel in open(UncompleteSpecLib):
    fullInput = pixel.split(',') #Split the text file on commas
    #string = fullInput[-1] #Save the last element
    #top = string.find('\n') #Find the newline character
    #fullInput[-1] = string[0:top] #Remove the newline character from the last element
    inPixel = fullInput[0] #Save the first element of list which should contain the polygon name that will be used to match
    #parse out the pixel's polygon name by finding the characters in front of the name and at end
    top = inPixel.find(")")
    bottom = (inPixel.find('=')+1)
    inPixel = inPixel[bottom:top]#Save just the pixels' polygon name

    j = 0 #set counter
    number = range(0,i)#create a list with numbers ranging from 0 to the length of the polygon metadata (first file that was looped through)
    for j in number: #Loop through polygon metadata
        if inPixel == meta[j][1]: #Check to see if the pixel belongs to the metadata polygon, if it does
            outputFile.write(fullInput[0] + ',')
            outputFile.write(', '.join(meta[j]))#Write to file the line of metadata ##+ fullInput[1:5]
            outputFile.write(',' + resolution +',' + date + ',' + sensor + ',' + flightline) #add additional data
            outputFile.write('\n')
            break
        j = j +1
outputFile.close()
#########################END#####################################################
    
