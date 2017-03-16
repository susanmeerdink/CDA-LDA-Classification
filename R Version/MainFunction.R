#####################################################################################################################
# Susan Meerdink
# CDA/LDA Classification
# 6/1/2016
# Main function, Main body, entry point
#####################################################################################################################
#### INPUTS ####
main_path <- 'F:\\Image-To-Image-Registration\\AVIRIS\\' #Set directory that holds all flightlines
shapefilePath <- 'C:\\Users\\Susan\\Documents\\GitHub\\CDA-LDA-Classification\\R Version\\Data Inputs\\2016_08_15_reference_polygons_withMeta_WGS84.shp'
#fl_list = ['FL02','FL03','FL04','FL05','FL06','FL07','FL08','FL09','FL10','FL11'] #Create the list of folders
fl_list <- {'FL03'} #Create the list of folders
#####################################################################################################################
#### Dependencies ####
library(rgdal)
library(raster)
require(maptools) 

#####################################################################################################################
#### Setting up Inputs
polygon <- shapefile(shapefilePath)#Load in Shapefile with polygons of plant species
proj4string(polygon)

#####################################################################################################################
#### Managing and Preparing Spectral Libraries ####

for (single_flightline in fl_list){ #Loop through flightlines
  print(paste('Starting with ',single_flightline))
  setwd(paste0(main_path,single_flightline,'\\')) #Set working directory
  image_list <- list.files(pattern = "*_ResizePlusBorder.dat") #Find all images with this suffix
  
  for (single_image in image_list){ #loop through files in Flightline
    single_path <- paste0(main_path,single_flightline,'\\',single_image) #Set path to a individual image
    raster <- brick(readGDAL(single_path))## Read in ENVI file and converts to raster brick, brick loads it into R rather than accessing the computers memory
    spectra <- data.frame(extract(raster, polygon),y=as.factor(polygon@data$Dominant_C)) #
    
  }
}

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