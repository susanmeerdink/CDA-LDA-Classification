directory = 'C:\\Users\\grad\\Dropbox\\Presentations\\2016_04_AAG\\Results'
producers = 'Producers_Accuracy.csv'
users = 'Users_Accuracy.csv'
overall = 'Overall_Accuracy.csv'
kappa = 'Kappa_Accuracy.csv'

import glob
import os
os.chdir(directory)

listFiles = glob.glob('*_Classification.csv')

producersFile = open(producers,'w')
usersFile = open(users,'w')
overallFile = open(overall,'w')
kappaFile = open(kappa,'w')

for oneFile in listFiles:
    end = oneFile.index('_OUTPUT')
    name = oneFile[0:end]
    
    for line in open(oneFile,'r'):
        lineAll = line.split(',')
        lineAll = map(lambda lineAll: lineAll.strip(), lineAll)
        
        if lineAll[0] == 'Producers Accuracy':
            producersFile.write(name + ',' + ','.join(lineAll)+ '\n')
        if lineAll[0] == 'Users Accuracy':
            usersFile.write(name + ',' + ','.join(lineAll)+ '\n')
        if lineAll[0] == 'Overall Accuracy':
            overallFile.write(name + ',' + ','.join(lineAll)+ '\n')
        if lineAll[0] == 'Kappa':
            kappaFile.write(name + ',' + ','.join(lineAll)+ '\n')
            
    
producersFile.close()
usersFile.close()
overallFile.close()
kappaFile.close()
