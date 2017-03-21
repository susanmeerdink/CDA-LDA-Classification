# Susan Meerdink
# last updated 3/18/2016
# This code reads in ROI formated csv (formatted using pullingOutROIInfo.py) and calculates
# classification results including producers, users, overal, and kappa accuracy
# results are output to a csv

#############INPUTS######################################
#inROI: should be your formatted ROI information with classification code results(formatted using pullingOutROIInfo.py)
#lookup: should be separate text file that contains in the first column classification code (0,1,2,..) and
#second column the class name (MAGF, ARGES, CEME,....) This column MUST match the metadata in the Dominant column inROI
#The classification code MUST be sequential, but it does not matter if it starts with 0 or 1
#inROI = 'R:\\users\\susan.meerdink\\Dropbox\\AAG_2016_Research\\ROIs\\20131125_AVIRIS&MASTER_same_OUTPUT.csv'
#inROI = 'R:\\users\\susan.meerdink\\Dropbox\\AAG_2016_Research\\ROIs\\20140829_AVIRIS&MASTER_diff_OUTPUT.csv'
#inROI = 'R:\\users\\susan.meerdink\\Dropbox\\AAG_2016_Research\\ROIs\\20140829_AVIRIS_same_OUTPUT.csv'
inROI = 'I:\\Classification-Products\\FL03\\3 - Classification Results\\f130411_class_results_OUTPUT_filled.csv'
lookup ='I:\\Classification-Products\\FL03\\3 - Classification Results\\Lookup_Table.csv'

############ENDINPUTS####################################

import numpy as np #Import numpy to work with arrays

#Open Files#
#output: will be the file that contains the classification results and accuracy as a .csv
end = inROI.index('.csv')
output = inROI[0:end] + '_Classification.csv'
outputFile = open(output,'w') #open the output file for writing

#Declare variables for lookup file
lookupNum = [] #VAriable that will hold the classification code (0, 1, 2, 3, ....)
lookupStr = [] #Variable that will hold the classification name (MAGF, AGRES, CEME, BAPI, ....)

#Read in Lookup File#
for line in open(lookup):
    lookupAll = line.split(',') #Split the line at commas into an array
    lookupAll = map(lambda lookupAll: lookupAll.strip(), lookupAll)
    lookupNum.append(int(lookupAll[0])) #Put the first column into the classification code variable
    lookupStr.append(lookupAll[1]) #Put the second column into the classification name variable

#Declare variables for classification results array
classArray = np.zeros((len(lookupNum),len(lookupNum))) #Create an array of zeros that is the number of class for row and columns
count = 0 #Used to tell if it beginning of file and to skip header
                      
#Read through input File (or formtted ROI file)#
for line in open(inROI):
    if count == 0: #If it is the header find column indexes
        lineAll = line.split(',') #Split the line at commas into an array
        lineAll = map(lambda lineAll: lineAll.strip(), lineAll)
        columnActual = lineAll.index('Dominant') #find the index of column labeled dominant which hold actual class
        columnModel = lineAll.index('Classification') #find index of column labeled with classification code                                   
    else: #If not read in the line
        lineAll = line.split(',') #Split the line at commas into an array
        lineAll = map(lambda lineAll: lineAll.strip(), lineAll)
        indexModel = int(lineAll[columnModel]) #Grab model classification code and thus index
        try:
            if indexModel != 0:   
                indexModel = indexModel -1  #to assign value to correct index
                indexActual = lookupStr.index(lineAll[columnActual]) #Find Actual code and thus index
                #indexActual = lookupStr.index(lineAll[columnActual].upper()) #Find Actual code and thus index - Uncomment this line for species codes

                #Assign Value#
                classArray[indexModel,indexActual] =  classArray[indexModel,indexActual] + 1 #Add one to this position
        except:
            continue

    count = count + 1 #Add value onto the count                      

#Summing up Values#
rowTotals = np.sum(classArray,axis = 1) #sum up each row
columnTotals = np.sum(classArray,axis = 0) #sum up each column
pixelTotal = np.sum(rowTotals) #Get total number of pixels

#Variables for Classification User's and Producer's Accuracy#
count = 0 #Set up value to loop through
producer = np.zeros((len(lookupNum))) #Create variable that will hold producer's accuracy
user = np.zeros((len(lookupNum))) #Create variable that will hold producer's accuracy
totalCorrect = 0 #Variable that will sum up total number of correctly classified pixels                      

#Classification User's and Producer's Accuracy#                      
while count < classArray.shape[0]: #While we haven't reached end of classification array
    classCorrect = classArray[count,count] #get the total number of pixels for class that have been classified correctly
    if classCorrect == 0:
        producer[count] = 0
        user[count] = 0
    else: 
        #Divide the correct number of pixels by the total number of actual pixels present                      
        producer[count] = classCorrect/columnTotals[count]
        #Divide the correct number of pixels by the total number of actual pixels present                      
        user[count] = classCorrect/rowTotals[count]

    totalCorrect = totalCorrect + classCorrect #add current correctly classified pixels to total correct pixels                      
    count = count + 1 #increase the counter
    
#Overall Accuracy#
overall = totalCorrect/pixelTotal #divide the total number of correctly classified pixels by the total number of pixels                      

#Kappa#
#topkappa: Multiply correctly classified pixels with total number of pixels subtract row totals
topKappa = (totalCorrect*pixelTotal) - np.sum(rowTotals, axis = 0)
#bottomkappa: Multiply total number of pixels with total number of pixels subtract row totals
bottomKappa = (pixelTotal*pixelTotal) - np.sum(rowTotals,axis = 0)
kappa = topKappa/bottomKappa #Get Kappa statistic by dividing top and bottom

#Write Results#
outputFile.write(',' + ','.join(lookupStr) + ',Row Total\n') #Write first line with header
#Loop through classification array
loopRange = range(0,classArray.shape[0])#Set range of values to loop through
for value in loopRange:
    outputFile.write(lookupStr[value]+',') #Write classification name
    outputFile.write(','.join(map(str,classArray[value,:]))) #Write out classification array
    outputFile.write(','+ str(rowTotals[value]) + '\n') #Finish row writing with row total

outputFile.write('Column Total,' + ','.join(map(str,columnTotals)) + ',' + str(pixelTotal) +'\n') #Write out column totals
outputFile.write('\n') #Add another blank row for visual
outputFile.write('Producers Accuracy,' + ','.join(map(str,producer)) + '\n') #Write out producer's accuracy
outputFile.write('Users Accuracy,' + ','.join(map(str,user)) + '\n') #Write out user's accuracy
outputFile.write('Number of Correctly Classified Pixels,' + str(totalCorrect) + '\n') #Write total number of correct pixels
outputFile.write('Overall Accuracy,' + str(overall) + '\n') #Write overall accuracy
outputFile.write('Kappa,' + str(kappa) + '\n') #Write out kappa statistic                                           

print('Finished Classification Accuracy Results found in:')
print(output)
#Close Files#
outputFile.close()                     
############ENDFILE####################################
