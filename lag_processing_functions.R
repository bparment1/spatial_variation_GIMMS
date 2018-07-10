############### General use: Lag processing fuction ########## 
## Generating lag profiles and lag raster from raster time series.
## General usage.
##
##
## DATE CREATED: 01/29/2018
## DATE MODIFIED: 07/10/2018
## AUTHORS: Benoit Parmentier  
## Version: 1
## PROJECT: spatial variability landscape
## ISSUE: 
## TO DO:
##
## COMMIT: 
##
## Links to investigate:


#### This files contains the following functions:

#1) autocor_filter_fun
#2) moran_multiple_fun
#3) local_moran_multiple_fun
#4) calculate_moranI_profile
#5) generate_lag_data_time_fun

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
  #3) multiband: if TRUE generate output as multiband file for each file in the stack
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
  
  #if(multiband==TRUE){
  #  d
  #}
  
  
  ### Read in filters used to computer Moran's I
  
  r <- subset(list_param$r_stack,i)
  
  if(is.null(NA_flag_val)){
    NA_flag_val<- NAvalue(r)
  }
  
  #moran_list <- MoranLocal(r,list_filters[[1]])
  moran_list <- lapply(1:length(list_filters),FUN=function(i,x){MoranLocal(x,w=list_filters[[i]])},x=r)
  r_local_moran <- stack(moran_list)

  n_lags <- length(list_filters) #number of lags considered
  names(r_local_moran) <- paste0("lag_",1:n_lags)
    
  plot(r_local_moran,y=1,main="lag 1",col=matlab.like(255))
  #plot(r_local_moran,y=5,main="lag 5",col=matlab.like(255))
  #title("lag 1")
  #moran_v <-as.data.frame(unlist(moran_list))
  #names(moran_v)<-names(r)
  #raster_name <- paste0(names(r),out_suffix,file_format)
  raster_name <- names(r)
  data_type_str <- dataType(r) #find the dataType, this should be a future input param
  
  if(is.null(NA_flag_val)){
      NA_flag_val <- NAvalue(r)
  }
    
  if(multiband==TRUE){
    if(out_suffix==""){
      out_suffix_s <- paste0("_lag","_",1,"_",n_lags,sep="")
    }else{
      out_suffix_s <- paste0("_lag","_",1,"_",n_lags,"_",out_suffix,sep="")
    }
    
    raster_name_tmp <- paste(raster_name,out_suffix_s,file_format,sep="")
    bylayer_val <- FALSE #don't write out separate layer files for each "band"
    out_raster_name <- paste(raster_name,"_",out_suffix_s,file_format,sep="") #files created
    
  }
  if(multiband==FALSE){
    bylayer_val <- TRUE #write out separate layer files for each "band"
    #### Attached to raster name output
    if(out_suffix==""){
      out_suffix_s <- names(r_local_moran)
    }else{
      out_suffix_s <- paste(names(r_local_moran),"_",out_suffix,sep="")
    }
    raster_name_tmp <- paste(raster_name,file_format,sep="") #for input in writeRaster
    out_raster_name <- paste(raster_name,"_",out_suffix_s,file_format,sep="") #files created
  }
  
  
  #bylayer_val <- T
  writeRaster(r_local_moran,
              filename=raster_name_tmp,
              bylayer=bylayer_val,
              suffix=out_suffix_s,
              overwrite=TRUE,
              NAflag=NA_flag_val,
              datatype=data_type_str,
              options=c("COMPRESS=LZW"))
  
  #out_raster_name <- 
  #### Prepare return object
  local_moran_obj <- list(r_local_moran,out_raster_name)
  names(local_moran_obj) <- c("r_local_moran","out_raster_name")
  
  return(local_moran_obj)
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

generate_lag_data_time_fun <- function(tile_index,grid_filename,r,multiband=T,file_format=".tif",num_cores=1,out_dir=NULL,out_suffix=NULL){
  
  if(!inherits(r,"Raster")){
    r <- stack(r)
  }
  
  #grid_filename <- out_tiles_filename
  tile_grid <- st_read(grid_filename)
  
  ### Plot tiles being processed:
  
  #tile_grid_selected <- as(tile_grid[tile_index,],"Spatial")
  tile_grid_selected <- tile_grid[tile_index,]
  tile_selected_spdf <- as(tile_grid_selected, "Spatial")
  
  plot(r,y=1)
  plot(tile_grid$geometry,add=T)
  plot(tile_grid_selected,add=T,col=NA,border="red")
  
  text(gCentroid(tile_selected_spdf),paste0("tile ",tile_index, " processed"))
  
  #text(st_centroid(tile_grid_selected),paste0("tile ",tile_index, " processed"))
  r_tile <- crop(r,tile_selected_spdf)
  
  #writeRaster(r_tile,)
  out_suffix <- out_suffix_s
  raster_name <- paste0("r_tile_",tile_index,"_",out_suffix_s,file_format)
  
  data_type_str <- dataType(r_tile) #find the dataTyre
  NA_flag_val_tmp <- NAvalue(r_tile)
  writeRaster(r_tile,
              filename=file.path(out_dir,raster_name),
              #bylayer=F,
              #suffix=paste(names(r),"_",out_suffix,sep=""),
              #format=format_raster,
              #suffix=paste(names(r)),
              overwrite=TRUE,
              NAflag=NA_flag_val_tmp,
              datatype=data_type_str,
              options=c("COMPRESS=LZW"))
  
  #generate filters for 10 lags: quick solution
  browser()
  
  list_filters<-lapply(1:10,
                       FUN=autocor_filter_fun,
                       f_type="queen") #generate 10 filters
  #moran_list <- lapply(list_filters,FUN=Moran,x=r)
  
  #r_stack <- r_tile
  
  list_param_moran <- list(list_filters=list_filters,
                           r_stack=r_tile,
                           multiband=multiband,
                           out_suffix=NULL,
                           out_dir=NULL)
  
  #moran_r <-moran_multiple_fun(1,list_param=list_param_moran)
  nlayers(r_tile) 
  
  strsplit(names(r),split="[.]")
  
  #debug(local_moran_multiple_fun)
  r_test <- local_moran_multiple_fun(1,list_param=list_param_moran)
  
  ### Change function to have option for multibands to reduce number of files outputs.
  
  #12:19
  local_moran_I_list <-mclapply(1:nlayers(r_tile), 
                                list_param=list_param_moran, 
                                FUN=local_moran_multiple_fun,
                                mc.preschedule=FALSE,
                                mc.cores = num_cores) #This is the end bracket from mclapply(...) statement
  
  #r_local_moran_stack <- stack(unlist(local_moran_I_list))
  
  ### Add extraction of profiles in sample:
  
  return(r_local_moran_stack)
}

######################### END OF SCRIPT ##############################

#n_layer <- length(modis_subset_layer_Day)

#if(n_layer==1){
#  r <- readGDAL(modis_subset_layer_Day) 
#  r  <-raster(r)
#}
#if(n_layer>1){
#  list_r <- lapply(modis_subset_layer_Day,function(x){raster(readGDAL(x))})
#  #r <- readGDAL(modis_subset_layer_Day) 
#  r <- stack(list_r)
#  get_names_layers <- function(x){val_extracted <- unlist(strsplit(x,":")); val_extracted[length(val_extracted)]}
#  names_layers <- unlist(lapply(modis_subset_layer_Day,get_names_layers))
#  names(r) <- names_layers
#}

#if(!is.null(scaling_factors)){ #if scaling factor exists, scale values...(not applied for QC flags!!!)
#  r <- scaling_factors[1]*r + scaling_factors[2]
#}
#Finish this part...write out
#names_hdf <- as.character(unlist(strsplit(x=basename(hdf_filename), split="[.]")))

#char_nb<-length(names_hdf)-2
#names_hdf <- names_hdf[1:char_nb]
#names_hdf <- paste(names_hdf,collapse="_") #this is the name of the hdf file with "." replaced by "_"
#raster_name <- paste(names_hdf,"_",out_suffix,file_format,sep="")
#out_dir_str <-  dirname(hdf)
#set output dir from input above
#if(n_layer==1){
#  raster_name_tmp <- raster_name
#  if(file_format==".tif"){
#    
#    writeRaster(r, 
#                NAflag=NA_flag_val,
#                filename=file.path(out_dir_s,raster_name_tmp),
#                bylayer=TRUE,
#                bandorder="BSQ",
#                overwrite=TRUE,
#                #datatype=data_type_str, #this should be a future option for reduced size!!!
#                options=c("COMPRESS=LZW")) #compress by default
#  }else{
#    
#    writeRaster(r, 
#                NAflag=NA_flag_val,
#                filename=file.path(out_dir_s,raster_name_tmp),
#                bylayer=TRUE,
#                bandorder="BSQ",
#                overwrite=TRUE)   
#  }
#}

#### Now deal with multiband e.g. MOD09 reflectance
#if(n_layer>1){
#  
#  #Write out as brick
#  data_type_str <- dataType(r) #find the dataType, this should be a future input param
#  if(is.null(NA_flag_val)){
#    NA_flag_val <- NAvalue(r)
#  }
#  
#  if(multiband==TRUE){
#    raster_name_tmp <- paste(names_hdf,"_",product_type,file_format,sep="")
#    bylayer_val <- FALSE #don't write out separate layer files for each "band"
#  }
#  if(multiband==FALSE){
#    raster_name_tmp <- raster_name
#    bylayer_val <- TRUE #write out separate layer files for each "band"
#  }
#  
#  if(file_format==".tif"){
#    writeRaster(r,
#                filename=file.path(out_dir_s,raster_name_tmp),
#                bylayer=bylayer_val,
#                #suffix=paste(names(r),"_",out_suffix,sep=""),
#                #format=format_raster,
#                suffix=paste(names(r)),
#                overwrite=TRUE,
#                NAflag=NA_flag_val,
#                datatype=data_type_str,
#                options=c("COMPRESS=LZW"))
#    
#  }else{
#    #Don't use compression option if not tif
#    writeRaster(r,
#                filename=file.path(out_dir_s,raster_name_tmp),
#                bylayer=multiband,
#                #suffix=paste(names(r),"_",out_suffix,sep=""),
#                #format=format_raster,
#                suffix=paste(names(r)),
#                overwrite=TRUE,
#                NAflag=NA_flag_val,
#                datatype=data_type_str)
#    
#  }
#  
#}