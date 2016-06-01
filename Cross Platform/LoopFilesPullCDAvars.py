directory = 'E:\meerdink\Dropbox\AAG_2016_Research\Spectral Libraries'
am = 'AVIRIS&MASTER_CDAvars.csv'
a = 'AVIRIS_CDAvars.csv'

import fnmatch
import os
os.chdir(directory)

#listFiles = glob.glob('*_CDAvars_CDAcoeffs.csv',recursive = True)
listFiles = []
for root, dirnames, filenames in os.walk(directory):
    for filename in fnmatch.filter(filenames, '*_CDAvars_CDAcoeffs.csv'):
        listFiles.append(os.path.join(root, filename))

amFile = open(am,'w')
aFile = open(a,'w')

for oneFile in listFiles:
    try:
        end = oneFile.index('_sort')
    except:
        end = oneFile.index('_CDAvars')
    start = oneFile.index('\\2')+1
    name = oneFile[start:end]
    lineTotal = []
    print('Reading in ' + name)
    
    for line in open(oneFile,'r'):
        lineAll = line.split(',')
        lineAll = map(lambda lineAll: lineAll.strip(), lineAll)
        lineTotal.append(lineAll[0])
    
    if 'AVIRIS&MASTER' in name: 
        amFile.write(name + ',' + ','.join(lineTotal)+ '\n')
    else:
        aFile.write(name + ',' + ','.join(lineTotal)+ '\n')


amFile.close()
aFile.close()

