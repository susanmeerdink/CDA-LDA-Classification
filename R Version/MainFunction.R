#####################################################################################################################
# Susan Meerdink
# CDA/LDA Classification
# 6/1/2016
# Main function, Main body, entry point
#####################################################################################################################
#### INPUTS ####
directory <- 'R:\\users\\susan.meerdink\\GitHub\\CDA-LDA-Classification\\R Version\\'
specLibLoc <- 'fl111107 r08&r09&r10&r11&r12&r14 spectral library.sli' #name of the ENVI spectral library containing spectral pulled from polygons in imagery 
metadataLoc <- 'fl111107 r08&r09&r10&r11&r12&r14 spectral library.csv' #Name of metadata associated with spectral library
#####################################################################################################################
#### Dependencies ####
##Loading R functions from RStoolbox in Github
RStoolboxFiles <- list.files(path = "R:\\users\\susan.meerdink\\GitHub\\RStoolbox\\R", pattern = '*.R',full.names = 'TRUE')
for (f in RStoolboxFiles){source(f)}

#####################################################################################################################
#### Managing and Preparing Spectral Libraries ####

## Using RStoolbox to read in ENVI spectral library
## This is a forked version from GitHub with the number of elements allowed increased.
## This function reads in an ENVI spectral library into a readable/workable format
specLib <- readSLI(path = paste(directory,specLibLoc,sep="")) 
print('Done reading in Spectral Library')

##read in the metadata for spectral library
metadata <- read.csv(paste(directory,metadataLoc,sep=""))
print('Done reading in metadata for Spectral Library')

##Strip redundant spectra from library
##In the process of resampling to 18 m using nearest, sometimes spectra are duplicated and should be removed.
##This function is originally based off Phil Dennison, Keely Roth, and Kenneth Dudley's code strip_redundant_spectra_v6.pro
specLibNoDup <- remove_duplicate_spectra(specLib,metadata)

#####################################################################################################################
#### Developing CDA Coefficients ####

#####################################################################################################################
#### LDA Classification ####

#####################################################################################################################
#### Accuracy Assessment ####