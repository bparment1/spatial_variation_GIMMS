##############################  Raster Processing Workflow  #######################################
#######################  Spatial utility: raster processing functions  ##############################
#This script contains utility functions to process raster data.
#
##AUTHOR: Benoit Parmentier 
#CREATED ON: 08/16/2018  
#MODIFIED ON: 08/22/2018            
#Version: 1
#PROJECT: General utilit     
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

## Get specific data type
#https://www.gdal.org/frmt_gtiff.html

generate_raster_dataType_table <- function(){
  
  ### Generate a table later on
  #Datatype definition	minimum possible value	maximum possible value
  
  vals <- c("LOG1S",NA,	FALSE,TRUE, 
            "INT1S",NA,	-127,	127,
            "INT1U",NA ,0, 255,
            "INT2S","Int16",	"-32,767","32,767",
            "INT2U",NA,	0,	"65,534",
            "INT4S","int32",	"-2,147,483,647",	"2,147,483,647",
            "INT4U",NA,	0,	"4,294,967,296",
            "FLT4S","Float32",	"-3.4e+38",	"3.4e+38",
            "FLT8S",NA,	"-1.7e+308",	"1.7e+308")
  
  
  dataType_table <- matrix(vals,nrow=9,ncol=4,byrow=T)
  
  dataType_table <-data.frame(dataType_table)
  
  names(dataType_table) <- c("r_type","gdal_type","min","max")
  ### bug error, columns have become factor: changed this here
  dataType_table <- data.frame(lapply(dataType_table, as.character), stringsAsFactors=FALSE)
  
  return(dataType_table)
}


############################ End of script ###################
