#####################################################################################################################
## Susan Meerdink
## 6/13/2016
## Remove duplicate spectra
## In the process of resampling to 18 m using nearest, sometimes spectra are duplicated and should be removed.
## This function is originally based off Phil Dennison, Keely Roth, and Kenneth Dudley's code strip_redundant_spectra_v6.pro
#####################################################################################################################
remove_duplicate_spectra <- function(input_spectral_lib,input_metadata){
  
  ##Misc Variables 
  tracker <-  0 #keeps track of how many duplicates are present in library
  dupIndex <- array() #keeps track of indices of duplicate spectra
  
  ##Find the sum of all reflectance values for each item
  sumLib <- colSums(input_spectral_lib)
  
  ##Loop through the library
  for(x1 in length(input_spectral_lib)){
    sum1 <- sumLib[x1]
    for(x2 in length(input_spectral_lib)){ ##Compare the above sum to the rest of the library
      if (x2 =! x1){ ##If the values does not equal the same index being compare continue processing
        sum2 <- sumLib[x2]
        if (sum1 == sum2){ ##If the sums are equal to each other, check to see if spectra are the same
          if (identical(input_spectral_lib[,x1],input_spectral_lib[,x2])){
            tracker <-  tracker + 1 #increase the count of duplicates
            append(dupIndex, x2) #add the index of duplicate to index list
          }
        }
      }
    }
  }
  
  ##Update Spectral Library
  output_spectral_library <- input_spectral_lib[,-c(dupIndex)] #remove duplicate indices
  
  ##Update Metadata
  output_metadata <- metadata[-c(dupIndex),]#remove duplicate indices
  
  ##Status/Results
  if (tracker == 0){
    print("There are no duplicate spectra in this library, quitting...")
  } else{
    print('There were ' + tracker + 'duplicate spectra in this library')
  }
  
  return(output_spectral_lib,output_metadata)
}