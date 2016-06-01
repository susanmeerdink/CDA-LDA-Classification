# Susan Meerdink
# last updated 12/8/2015
# ############INPUTS######################################
#MUST add semi colon (;) to end of file!
OrigROI = 'C:\Users\grad\Dropbox\Code\CDA_Work_With_Keely\Results\\AVIRIS_class_ROI_Output_formatted.csv'
MasterMetaData = 'C:\Users\grad\Dropbox\Code\CDA_Work_With_Keely\Results\\MountainPolyList.txt'
OutputROI = 'C:\Users\grad\Dropbox\Code\CDA_Work_With_Keely\Results\\AVIRIS_class_ROI_Output_Mountain.csv'
############ENDINPUTS####################################

import numpy as np
#Set Output file for ROIs
outputFile = open(OutputROI,'w')
#inputFile = open(OrigROI,'r')


valList = [] #Variable to hold metadata
#Loop through file that contains all the metadata fields for each polygon and store the metadata so that it can be assigned to the new library
for polygon in open(MasterMetaData):
    text = polygon.strip()#Split the line into array based on comma
    valList.append(text) #Add this line (polygon info) to the meta data variable 

i = 0
for line in open(OrigROI):
    if i ==0:  
        #Create headers for the output file
        outputFile.write(line)
    else:
        inLine = line.split(',')
        #print inLine
        if inLine[0] in valList:
            #print inLine[0]
            outputFile.write(','.join(inLine))
    i = i +1
    
outputFile.close()
#inputFile.close()
    

#########################END#####################################################
    
