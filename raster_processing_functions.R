##############################  Raster Processing Workflow  #######################################
#######################  Spatial utility: raster processing functions  ##############################
#This script contains utility functions to process raster data.
#
##AUTHOR: Benoit Parmentier 
#CREATED ON: 08/16/2018  
#MODIFIED ON: 08/24/2018            
#Version: 1
#PROJECT: General utility     
#COMMENTS:
#
#COMMIT: 
#
#TODO:
#1) 

#################################################################################################

### Loading R library and packages        
#library used in the workflow production:
library(gtools)                              # loading some useful tools 
library(mgcv)                                # GAM package by Simon Wood
library(sp)                                  # Spatial pacakge with class definition by Bivand et al.
library(spdep)                               # Spatial pacakge with methods and spatial stat. by Bivand et al.
library(rgdal)                               # GDAL wrapper for R, spatial utilities
#library(gstat)                               # Kriging and co-kriging by Pebesma et al., not on NEX
library(fields)                              # NCAR Spatial Interpolation methods such as kriging, splines
library(raster)                              # Hijmans et al. package for raster processing
library(gdata)                               # various tools with xls reading, cbindX
library(rasterVis)                           # Raster plotting functions
library(parallel)                            # Parallelization of processes with multiple cores
library(maptools)                            # Tools and functions for sp and other spatial objects e.g. spCbind
library(maps)                                # Tools and data for spatial/geographic objects
library(reshape)                             # Change shape of object, summarize results 
library(plotrix)                             # Additional plotting functions
library(plyr)                                # Various tools including rbind.fill
library(spgwr)                               # GWR method, not on NEX
#library(automap)                             # Kriging automatic fitting of variogram using gstat, not on NEX
library(rgeos)                               # Geometric, topologic library of functions
library(RPostgreSQL)                         # Interface R and Postgres, not used in this script
library(gridExtra)
#Additional libraries not used in workflow
library(pgirmess)                            # Krusall Wallis test with mulitple options, Kruskalmc {pgirmess}  
library(colorRamps)
library(zoo)
library(xts)

#### FUNCTIONS

split_multiband <- function(in_file,out_suffix_str,out_dir){
  
  ### Can also use gdal merge to split?
  
  r_in <- brick(in_file)
  
  #if data is multiband
  n_layers <- nlayers(r_in)
  #bylayer=F
  suffix_str <- paste0(1:n_layers,"_",out_suffix_str)
  
  bylayer_val <- T
  #raster_name_tmp <- "test.tif"
  file_format <- extension(in_file)
  #out_suffix
  
  #Split and get a raster stack
  r_out <- writeRaster(r_in,
                       filename=file.path(out_dir,basename(in_file)),
                       bylayer=bylayer_val,
                       suffix=suffix_str,
                       overwrite=TRUE,
                       #NAflag=NA_flag_val,
                       #datatype=data_type_str,
                       options=c("COMPRESS=LZW"))
  
  #
  lf_out <- file.path(out_dir,
                      paste0(names(r_out),file_format))
  
  return(lf_out)
} 

## 

generate_raster_dataType_table <- function(){
  #Goal: this function generate a table (data.frame) with data types
  # and valid value range used in the raster package R. The corresponding
  # data type in the GDAL library is provided to allow matching when using
  # GDAL commands.
  
  # Note that we are using the specific data types for tif.
  # The following links provide more information:
  #https://www.gdal.org/frmt_gtiff.html
  #urrently band types of Byte, UInt16, Int16, UInt32, Int32, Float32, 
  #Float64, CInt16, CInt32, CFloat32 and CFloat64 are supported for reading and writing.
  
  ######### Start scripts ################
  
  vals <- c("LOG1S",NA,	FALSE,TRUE, 
            "INT1S","Byte",	-127,	127,
            "INT1U",NA,0, 255,
            "INT2S","Int16",	"-32,767","32,767",
            "INT2U","UInt16",	0,	"65,534",
            "INT4S","int32",	"-2,147,483,647",	"2,147,483,647",
            "INT4U","UInt32",	0,	"4,294,967,296",
            "FLT4S","Float32",	"-3.4e+38",	"3.4e+38",
            "FLT8S","Float64",	"-1.7e+308",	"1.7e+308")
  
  dataType_table <- matrix(vals,nrow=9,ncol=4,byrow=T)
  
  dataType_table <-data.frame(dataType_table)
  
  names(dataType_table) <- c("r_type","gdal_type","min","max")
  ### bug error, columns have become factor: changed this here
  dataType_table <- data.frame(lapply(dataType_table, as.character), stringsAsFactors=FALSE)
  
  #class(dataType_table$gdal_type)
  
  return(dataType_table)
}



grid_sampling_raster <- function(x_sampling,y_sampling,r){
  grid_sf<- st_make_grid(r,n=c(x_sampling,y_sampling))
  grid_samples <- st_centroid(grid_sf)
  #plot(grid_sf,add=T)
  #plot(r,y=1)
  #plot(grid_samples,add=T)
  #class(grid)
  return(grid_samples)
}

############################ End of script ###################
