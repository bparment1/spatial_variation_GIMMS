############### SESYNC Research Support: Hurricane Management ########## 
## Help in processing modis data for the workshop group.
##
##
## DATE CREATED: 12/13/2017
## DATE MODIFIED: 12/14/2017
## AUTHORS: Benoit Parmentier  
## Version: 1
## PROJECT: spatial variability landscape
## ISSUE: 
## TO DO:
##
## COMMIT: 
##
## Links to investigate:

###################################################
#

###### Library used

library(gtools)                              # loading some useful tools 
library(sp)                                  # Spatial pacakge with class definition by Bivand et al.
library(spdep)                               # Spatial pacakge with methods and spatial stat. by Bivand et al.
library(raster)                              # raster functions and spatial utilities
library(rgdal)                               # GDAL wrapper for R, spatial utilities
library(gdata)                               # various tools with xls reading, cbindX
library(parallel)                            # Parallelization of processes with multiple cores
library(maptools)                            # Tools and functions for sp and other spatial objects e.g. spCbind
library(maps)                                # Tools and data for spatial/geographic objects
library(plyr)                                # Various tools including rbind.fill
library(dplyr)                               # data wrangling
library(rgeos)                               # Geometric, topologic library of functions
library(gridExtra)                           # Combining lattice plots
library(colorRamps)                          # Palette/color ramps for symbology
library(ggplot2)                             # plotting functionality

###### Functions used in this script and sourced from other files

create_dir_fun <- function(outDir,out_suffix=NULL){
  #if out_suffix is not null then append out_suffix string
  if(!is.null(out_suffix)){
    out_name <- paste("output_",out_suffix,sep="")
    outDir <- file.path(outDir,out_name)
  }
  #create if does not exists
  if(!file.exists(outDir)){
    dir.create(outDir)
  }
  return(outDir)
}

#Used to load RData object saved within the functions produced.
load_obj <- function(f){
  env <- new.env()
  nm <- load(f, env)[1]
  env[[nm]]
}

###### Functions used in this script

#Benoit setup
script_path <- "/nfs/bparmentier-data/Data/projects/spatial_variation_GIMMS/scripts"

#raster_processing_functions <- "raster_processing_functions_12142017.R" #Functions used to mosaic predicted tiles
#source(file.path(script_path,raster_processing_functions)) #source all functions used in this script 

############################################################################
#####  Parameters and argument set up ###########

#ARGS 1
#Set up teamhurricane-data
in_dir <- "/nfs/bparmentier-data/Data/projects/spatial_variation_GIMMS/data"

out_dir <- "/nfs/bparmentier-data/Data/projects/spatial_variation_GIMMS/outputs"

## Set up Benoit
#in_dir <- "/nfs/bparmentier-data/Data/projects/managing_hurricanes/data"
#out_dir <- "/nfs/bparmentier-data/Data/projects/managing_hurricanes/outputs"

hdf_file <-"test.hdf"

#NA_flag <- -999999
file_format <- ".tif"
scaling_factor <- 0.0001 #MODIFY THE SCALING FACTOR - FOR NORMALIZED DATA SHOULD BE 10,000 AT LEAST
#ARGS 7
create_out_dir_param=TRUE #create a new ouput dir if TRUE
#ARGS 8
out_suffix <-"GIMMS_processing_01242018" #output suffix for the files and ouptut folder #param 12
num_cores <- 2 # number of cores

################# START SCRIPT ###############################

######### PART 0: Set up the output dir ################

options(scipen=999)

if(is.null(out_dir)){
  out_dir <- in_dir #output will be created in the input dir
  
}
#out_dir <- in_dir #output will be created in the input dir

out_suffix_s <- out_suffix #can modify name of output suffix
if(create_out_dir_param==TRUE){
  out_dir <- create_dir_fun(out_dir,out_suffix)
  setwd(out_dir)
}else{
  setwd(out_dir) #use previoulsy defined directory
}

#######################################
### PART I READ AND PREPARE DATA #######
#set up the working directory
#Create output directory

## Download data here:
#https://ecocast.arc.nasa.gov/data/pub/gimms/3g.v1/
lf <- list.files(dir="https://ecocast.arc.nasa.gov/data/pub/gimms/3g.v1/")

lf_name <- download.file("https://ecocast.arc.nasa.gov/data/pub/gimms/3g.v1/00FILE-LIST.txt",
                         destfile = "00FILE-LIST.txt")

lf_df <- read.table("00FILE-LIST.txt",stringsAsFactors = F)

file1 <- download.file(lf_df[1,1],basename(lf_df[1,1]))
#https://ecocast.arc.nasa.gov/data/pub/gimms/
  
raster_file <- basename(lf_df[1,1])

GDALinfo_raster <- GDALinfo(raster_file,returnScaleOffset = F) #use GDAL info utility


#hdf_file <- file.path(in_dir,hdf_file) #associate file path to input
#GDALinfo_hdf <- GDALinfo(hdf_file,returnScaleOffset = F) #use GDAL info utility

#str(GDALinfo_hdf) ## Class GDLobj
raster_subdataset <- attributes(GDALinfo_raster)$subdsmdata
#print(modis_subdataset)

### Generate a data.frame from scientific datasets
#hdf_df <- (strsplit(modis_subdataset,":"))
#length(modis_subdataset)
#hdf_df <- as.data.frame(do.call(rbind,hdf_df),stringsAsFactors=F)
#hdf_df <- data.frame(lapply(hdf_df, as.character))
#hdf_df %>% 
#  mutate_all(as.character)

#names(hdf_df) <- c("subdataset_name","description","dir","product","var_name")
#Select automatically QC flag!!
#View(hdf_df)

#write.table(hdf_df,"hdf_subdataset.txt",sep=",")

#modis_subset_layer_Day <- paste("HDF4_EOS:EOS_GRID:",
#                                hdf_file,subdataset,sep="")

#NDVI variable
#modis_layer_str1 <- unlist(strsplit(modis_subdataset[1],"\""))[3] #Get day NDVI layer
#QC
#modis_layer_str2 <- unlist(strsplit(modis_subdataset[5],"\""))[3] #Get day VI QC layer

#subdataset <- modis_layer_str1
#modis_subset_layer_Day <- paste("HDF4_EOS:EOS_GRID:",hdf_file,subdataset,sep="")

#r <- readGDAL(modis_subset_layer_Day) #read specific dataset in hdf file and make SpatialGridDataFrame

#r  <-raster(r) #convert to raser object

#plot(r,main="NDVI ~250m")
#r #print properties
#res(r) #spatial resolution
#NAvalue(r) #find the NA flag value
#dataType(r) #find the dataTyre

# for more control, you can set dataType and/or compress the files
#data_type_str <- "FLT4S"
#NA_flag_val <- NAvalue(r2)

#writeRaster(r2,
#            file.path(out_dir,raster_name),
#            overwrite=TRUE,
#            NAflag=NA_flag_val,
#            datatype=data_type_str,
#            options=c("COMPRESS=LZW"))

#### Next steps to consider:
## Use the name from MODIS file because it contains information on tile location, date and product type
## Use QC index to screen for low value pixels

######################### END OF SCRIPT ##############################