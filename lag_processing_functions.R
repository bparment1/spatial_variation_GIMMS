############### General use: Lag processing fuction ########## 
## Generating lag profiles and lag raster from raster time series.
## General usage.
##
##
## DATE CREATED: 01/29/2018
## DATE MODIFIED: 06/28/2018
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
library(sf)

###### Functions used in this script and sourced from other files


### generate filter for Moran's I function in raster package
autocor_filter_fun <-function(no_lag=1,f_type="queen"){
  if(f_type=="queen"){
    no_rows <- 2*no_lag +1
    border_row <-rep(1,no_rows)
    other_row <- c(1,rep(0,no_rows-2),1)
    other_rows <- rep(other_row,no_rows-2)
    mat_data<- c(border_row,other_rows,border_row)
    autocor_filter<-matrix(mat_data,nrow=no_rows)
  }
  #if(f_type=="rook){} #add later
  return(autocor_filter)
}


#MODIFY: calculate for multiple dates and create averages...
#Now run Moran's I for raster image given a list of  filters for different lags and raster stack
moran_multiple_fun<-function(i,list_param){
  #Parameters:
  #list_filters: list of filters with different lags in the image
  #r_stack: stack of raster image, only the selected layer is used...
  list_filters <-list_param$list_filters
  r <- subset(list_param$r_stack,i)
  moran_list <- lapply(list_filters,FUN=Moran,x=r)
  moran_v <-as.data.frame(unlist(moran_list))
  names(moran_v)<-names(r)
  return(moran_v)
}

local_moran_multiple_fun<-function(i,list_param){
  #
  #INPUTS:
  # i: select raster file to process
  #1) list_filters: list of filters with different lags in the image
  #2) r_stack: stack of raster image, only the selected layer is used...
  #3) NA_flag_val: if NULL, it is derived from input stack
  #4) file_format: default is *.tif (if NULL)
  #5) out_suffix: if NULL default, no suffix added.
  #6) out_dir: if NULL default, use current dir
  
  #OUTPUTS:
  #
  
  ###### Begin #######
  
  
  ### First, read in all parameters:
  
  list_filters <- list_param$list_filters
  list_r_stack <- list_param$r_stack
  NA_flag_val <- list_param$NA_flag_val
  file_format <- list_param$file_format
  out_suffix <- list_param$out_suffix
  out_dir <- list_param$out_dir
  
  #### Now check default values:
  
  if(is.null(out_suffix)){
    out_suffix <- ""
  }
  if(is.null(out_dir)){
    out_dir <- "."
  }
  
  ### Read in filters used to computer Moran's I
  
  r <- subset(list_param$r_stack,i)
  
  if(is.null(NA_flag_val)){
    NA_flag_val<- NAvalue(r)
  }
  
  #moran_list <- MoranLocal(r,list_filters[[1]])
  moran_list <- lapply(1:length(list_filters),FUN=function(i,x){MoranLocal(x,w=list_filters[[i]])},x=r)
  r_local_moran <- stack(moran_list)
  
  names(r_local_moran) <- paste0("lag_",1:length(list_filters))
  
  plot(r_local_moran,y=1,main="lag 1",col=matlab.like(255))
  #plot(r_local_moran,y=5,main="lag 5",col=matlab.like(255))
  #title("lag 1")
  #moran_v <-as.data.frame(unlist(moran_list))
  #names(moran_v)<-names(r)
  raster_name <- paste0(names(r),out_suffix,file_format)
  data_type_str <- dataType(r)
  
  
  #### Attached to raster name output
  if(out_suffix==""){
    out_suffix_s <- names(r_local_moran)
  }else{
    out_suffi_s <- paste(names(r_local_moran),"_",out_suffix,sep="")
  }
  
  writeRaster(r_local_moran,
              filename=raster_name,
              bylayer=T,
              suffix=out_suffix_s,
              overwrite=TRUE,
              NAflag=NA_flag_val,
              datatype=data_type_str,
              options=c("COMPRESS=LZW"))
  
  #### Prepare return object
  
  return(r_local_moran)
}

#Extract moran's I profile from list of images...the list may contain sublist!!! e.g. for diffeferent
#methods in interpolation
calculate_moranI_profile <- function(lf,nb_lag){
  list_filters<-lapply(1:nb_lag,FUN=autocor_filter_fun,f_type="queen") #generate lag 10 filters
  #moran_list <- lapply(list_filters,FUN=Moran,x=r)
  list_moran_df <- vector("list",length=length(lf))
  for (j in 1:length(lf)){
    r_stack <- stack(lf[[j]])
    list_param_moran <- list(list_filters=list_filters,r_stack=r_stack) #prepare parameters list for function
    #moran_r <-moran_multiple_fun(1,list_param=list_param_moran)
    nlayers(r_stack) 
    moran_I_df <-mclapply(1:nlayers(r_stack), list_param=list_param_moran, FUN=moran_multiple_fun,mc.preschedule=FALSE,mc.cores = 10) #This is the end bracket from mclapply(...) statement
    
    moran_df <- do.call(cbind,moran_I_df) #bind Moran's I value 10*nlayers data.frame
    moran_df$lag <-1:nrow(moran_df)
    
    list_moran_df[[j]] <- moran_df
  }
  names(list_moran_df) <- names(lf)
  return(list_moran_df)
}

######################### END OF SCRIPT ##############################

