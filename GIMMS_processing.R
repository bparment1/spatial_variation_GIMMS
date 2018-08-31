############### SESYNC Research Support: GIMMS Spatial Variation ########## 
## Processing GIMMS NDVI data to explore spatio-temporal variability.
## GIMMS covers the 1982-2015 time period at 0.08333 degree resolution (~8km).
## General processing function for climatelandfeedback
##
## DATE CREATED: 01/24/2018
## DATE MODIFIED: 08/31/2018
## AUTHORS: Benoit Parmentier  
## Version: 1
## PROJECT: spatial variability landscape
## ISSUE: 
## TO DO:
##
## COMMIT: generating table for data type and range for Raster and GDAL formats
##
##Data downloaded from:
#https://ecocast.arc.nasa.gov/data/pub/gimms/3g.v1/

###################################################
#

###### Library used

library(sp) # spatial/geographic objects and functions
library(rgdal) #GDAL/OGR binding for R with functionalities
library(spdep) #spatial analyses operations, functions etc.
library(gtools) # contains mixsort and other useful functions
library(maptools) # tools to manipulate spatial data
library(parallel) # parallel computation, part of base package no
library(rasterVis) # raster visualization operations
library(raster) # raster functionalities
library(forecast) #ARIMA forecasting
library(xts) #extension for time series object and analyses
library(zoo) # time series object and analysis
library(lubridate) # dates functionality
library(colorRamps) #contains matlab.like color palette
library(rgeos) #contains topological operations
library(sphet) #contains spreg, spatial regression modeling
library(BMS) #contains hex2bin and bin2hex, Bayesian methods
library(bitops) # function for bitwise operations
library(foreign) # import datasets from SAS, spss, stata and other sources
library(gdata) #read xls, dbf etc., not recently updated but useful
library(classInt) #methods to generate class limits
library(plyr) #data wrangling: various operations for splitting, combining data
#library(gstat) #spatial interpolation and kriging methods
library(readxl) #functionalities to read in excel type data
library(psych) #pca/eigenvector decomposition functionalities
library(snow)
library(sf)

###### Functions used in this script and sourced from other files

#Benoit setup
script_path <- "/nfs/bparmentier-data/Data/projects/spatial_variation_GIMMS/scripts"

gimms_processing_functions <- "GIMMS_processing_functions_07252018.R" #Functions used to mosaic predicted tiles
generate_tiles_functions <- "generate_spatial_tiles_functions_07102018.R" #Functions used to mosaic predicted tiles
lag_processing_functions <- "lag_processing_functions_07302018.R"
get_study_region_functions <- "get_study_region_data_07252018.R"
mosaicing_functions <- "weighted_mosaicing_functions_08302018.R"
raster_processing_functions <- "raster_processing_functions_08242018.R"
source(file.path(script_path,gimms_processing_functions)) #source all functions used in this script 
source(file.path(script_path,generate_tiles_functions))
source(file.path(script_path,lag_processing_functions))
source(file.path(script_path,get_study_region_functions))
source(file.path(script_path,mosaicing_functions))
source(file.path(script_path,raster_processing_functions))

#########cd ###################################################################
#####  Parameters and argument set up ########### 

#ARGS 1
in_dir <- "/nfs/bparmentier-data/Data/projects/spatial_variation_GIMMS/data"
#ARGS 2
out_dir <- "/nfs/bparmentier-data/Data/projects/spatial_variation_GIMMS/outputs"
#ARGS 3:
#NA_flag <- -999999
NA_flag_val <- NULL
#ARGS 4:
file_format <- ".tif"
#ARGS 5:
scaling_factor <- 0.0001 #MODIFY THE SCALING FACTOR - FOR NORMALIZED DATA SHOULD BE 10,000 AT LEAST
#ARGS 6
create_out_dir_param=TRUE #create a new ouput dir if TRUE
#ARGS 7
out_suffix <-"GIMMS_processing_08302018" #output suffix for the files and ouptut folder
#ARGS 8
num_cores <- 2 # number of cores
#ARGS 9
date_param <- "1982.01.01;1982.12.31" #start date, end date
#ARGS 10
GIMMS_product <- "3g.v1"

#ARGS 11
#we want 20x20
tile_ratio <- c(360/20,180/20) # in the order for x and y
#ARGS 12:generate overlapping tiles
tile_overlap <- c(0.25,0.25) # in the order for x and y
#ARGS 13
multiband <- TRUE
#ARGS 14
max_lag <- 10 #maximum lag to consider

#ARGS 15
processing_steps <- list(download=FALSE,
                         import=FALSE,
                         tiling=TRUE)

################# START SCRIPT ###############################

######### PART 0: Set up the output dir ################

options(scipen=999)

#set up the working directory
#Create output directory

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
### PART 1: DOWNLOAD DATA #######

#Will be a function
## Make this a function later on!!!
date_param <- unlist(strsplit(date_param,";"))

start_date <- date_param[1]
end_date <- date_param[2]

if(processing_steps$download==TRUE){
  
  
  #debug(gimms_product_download)
  #undebug(GIMMS_product_download)
  
  #raster_processing_functions <- "GIMMS_processing_functions_02052018.R" #Functions used to mosaic predicted tiles
  #source(file.path(script_path,raster_processing_functions)) #source all functions used in this script 
  
  gimms_download_obj <- GIMMS_product_download(GIMMS_product,
                                 start_date,
                                 end_date,
                                 out_dir=in_dir, #store nc4 in data folder
                                 out_suffix)  ##Functions used in the script
  f <- gimms_download_obj$downloaded_files
  
  
}else{
  f <- list.files(path=in_dir,
                  pattern="*.nc4",
                  full.names=T)
}


##############################
#### PART 2: Import dataset
####
##### Next import in Tif format the NCDF

if(processing_steps$import==TRUE){
  
  #f <- test$downloaded_files
  #f <- list.files(path=in_dir,pattern="*.nc4",full.names=T)
  
  #debug(import_gimms_nc4)
  
  lf_imported <- import_gimms_nc4(input_file=f[1],
                                  var_name="ndvi",
                                  NA_flag_val=NULL,
                                  file_format=file_format,
                                  out_suffix="",
                                  out_dir=in_dir)
  
  #input_file,var_name="ndvi",NA_flag_val,file_format,out_suffix="",out_dir="."
  
  lf_imported <- imported_files <- mclapply(f,
                                            FUN=import_gimms_nc4,
                                            var_name="ndvi",
                                            NA_flag_val=NULL,
                                            file_format=file_format,
                                            out_suffix="",
                                            out_dir=in_dir,
                                            mc.preschedule=FALSE,
                                            mc.cores = num_cores)
  
}else{
  f <- list.files(path=in_dir,
                  pattern=paste0("*",file_format),
                  full.names=T)
}

#### Next steps to consider:
## Use QC index to screen for low value pixels

#########################################
########### PART 3: Generate grid ###########

#This can be non-overlapping or overlapping

lf_gimms <- mixedsort(f)

ref_file <- lf_gimms[1]

##### Generate a grid/tile for processing:
## Must transformed to a function later on.

if(processing_steps$tiling==TRUE){
  
  r_ref <- raster(ref_file)
  
  plot(r_ref)
  
  ### store in the tile directory
  #undebug(generate_tiles_from_extent)
  
  #we want 20x20
  x_ratio <- tile_ratio[1]
  y_ratio <- tile_ratio[2]
  
  ### generate overlapping tiles
  x_overlap <- tile_overlap[1]
  y_overlap <- tile_overlap[2]
  
  #debug(generate_tiles_from_extent)
  tiles_obj <- generate_tiles_from_extent(r_ref,
                                          y_ratio=y_ratio,
                                          x_ratio=x_ratio,
                                          y_overlap=y_overlap,
                                          x_overlap=x_overlap,
                                          out_suffix=out_suffix,
                                          out_dir=NULL) #this does not work...
  
  df_tiles <- tiles_obj$df_tiles

  #### Transform a data.frame object to a sf object using coordinates x and y
  centroids_tiles <- st_as_sf(df_tiles,
                      coords = c('x_center', 'y_center'),
                      crs = proj4string(r_ref))
  plot(r_ref)
  text(df_tiles$x_center,df_tiles$y_center,df_tiles$ID,cex=0.5)
  centroids_tiles$ID
  
  out_centroids_filename <- file.path(out_dir,"tiles_centroids.shp")
  st_write(centroids_tiles,
           dsn=out_centroids_filename,
           delete_dsn =T)
  
  tiles_combined_spdf <- do.call(rbind,
                                 tiles_obj$list_tiles)
  
  plot(tiles_combined_spdf)
  tiles_combined_sf <- as(tiles_combined_spdf,"sf") #transforming spdf to sf object
  st_crs(tiles_combined_sf) <- proj4string(r_ref)
  
  tiles_combined_sf$ID <- 1:nrow(tiles_combined_sf)
  
  #tiles_all <- join(tiles_combined_sf,df_tiles)
  tiles_combined_sf <- merge(tiles_combined_sf,df_tiles,by="ID")
  #View(tiles_combined_sf)

  out_tiles_filename <- file.path(out_dir,"tiles_combined.shp")
  st_write(tiles_combined_sf,
           dsn=out_tiles_filename,
           delete_dsn = T)
  
  ############# This part should be a function: select only the tiles that are within this area:
  
  country_names <- c("Kazakhstan") #can have mulitple names as a vector here
  country_names <- c("Mongolia")
  
  #undebug(get_countries_outlines)
  reg_outline_spdf <- get_countries_outlines(country_names,
                                             out_dir=".",
                                             out_suffix=country_names)
  reg_outline_sf <- as(reg_outline_spdf,"sf")
  #plot(reg_outline_sf)
  selected_poly_ID <- st_intersects(reg_outline_sf,tiles_combined_sf)
  selected_poly_ID <- unlist(selected_poly_ID)
  
  tiles_reg_sf <- filter(tiles_combined_sf,ID%in%selected_poly_ID)
  
  ## Let's drop 41 since the overlap is low
  tiles_reg_sf <- filter(tiles_reg_sf,ID!=41)
  ### visualize selected polygons
  plot(r_ref)
  text(df_tiles$x_center,df_tiles$y_center,df_tiles$ID,cex=0.5)
  plot(reg_outline_sf$geometry,add=T,border="blue")
  plot(tiles_reg_sf$geometry,border="red",add=T)
  
  out_tiles_filename <- file.path(out_dir,"tiles_reg.shp")
  st_write(tiles_reg_sf,
           dsn=out_tiles_filename,
           delete_dsn = T)
  
  ##### end of selected file
  
}


#########################################
########### PART 4: Analysis of spatial variability ###########

#This part should be a separate script to run on cluster

### Create function for SLURM job:

#ARGS 11
#infile_list_tiles <- "list_tiles.txt"
#ARGS 12: 
#TILE INDEX this is from the array
# Get the the array id value from the environment variable passed from sbatch

#SLURM_ARRAY_TASK_ID <- Sys.getenv('SLURM_ARRAY_TASK_ID')
#tile_index <- Sys.getenv("SLURM_ARRAY_TASK_ID") #this is should be an integer from 1:n, n is the number of tiles
#slurm_arrayid <- Sys.getenv('SLURM_ARRAYID') #work with #SARRAY option in SLURM
#tile_index <- as.numeric(slurm_arrayid) # coerce the value to an integer
#tile_index <- 1  #for testing

#### Arguments for function
#tile_index
#grid_filename
#r_reg: raster file stack (can be one image too) or one image
#multiband
#out_dir
#out_suffix

tile_index <- 4 #This is tile grid 65
lf_gimms <- mixedsort(list.files(pattern=file_format,path=in_dir,full.names=T))

tiles_sf <- st_read(out_tiles_filename)
n_tiles <- nrow(tiles_sf)
#ref_file <- stack(lf_gimms[1]

#debug(generate_lag_data_time_fun)

test <- generate_lag_data_time_fun(tile_index=tile_index,
                                   grid_filename=out_tiles_filename,
                                   r=lf_gimms,
                                   max_lag=max_lag,
                                   multiband=multiband,
                                   file_format=file_format,
                                  num_cores=4,
                                   out_dir=NULL,#if null a output dir name with tile_+index nunber is created
                                   out_suffix=out_suffix)

#undebug(generate_lag_data_time_fun)
#test <- lapply(1:n_tiles,
#               FUN=generate_lag_data_time_fun,
#               grid_filename=out_tiles_filename,
#               r=lf_gimms,
#               max_lag=max_lag,
#               multiband=multiband,
#               file_format=file_format,
#              num_cores=4,
#               out_dir=NULL,
#               out_suffix=out_suffix)

test <- mclapply(1:n_tiles,
               FUN=generate_lag_data_time_fun,
               grid_filename=out_tiles_filename,
               r=lf_gimms,
               max_lag=max_lag,
               multiband=multiband,
               file_format=file_format,
               num_cores=4,
               out_dir=NULL,
               out_suffix=out_suffix)

#test[[1]][[1]]
#test[[1]][[1]]$out_raster_name

#> names(test[[1]][[1]])
#[1] "r_local_moran"   "out_raster_name"
#> (test[[1]][[1]]$out_raster_name)
#[1] "ndvi3g_geo_v1_NDVI_1982_01_01__lag_1_10.tif"
#> (test[[1]][[2]]$out_raster_name)
#[1] "ndvi3g_geo_v1_NDVI_1982_01_15__lag_1_10.tif"

length(test[[1]][[1]]) #24

length(test) #8
### Reorganize the files for mosaicing

test[[1]][1]
##### Now mosaic:

data_type_str <- dataType(test[[1]][[1]]$r_local_moran)


#########################################
########### PART 5: Mosaic outputs ###########

if(processing_steps$mosaicing==TRUE){
  
  ### PARAM to prepare for mosaicing
  
  #1) lf_mosaic
  #2) mosaic_method="unweighted", #default val
  #3) num_cores=1,
  #4) r_mask_raster_name<-NULL  #this gives the extent
  #5) python_bin=NULL
  #6) mosaic_python="/nobackupp6/aguzman4/climateLayers/sharedCode/gdal_merge_sum_noDataTest.py",
  #7) algorithm="R",
  #8) match_extent=TRUE,
  #9) df_points=NULL,
  #10) NA_flag_val= HULL #-9999,
  #11) file_format=".tif" #default val
  #12) out_suffix=NULL,
  #13) out_dir=NULL,
  #14) tmp_files=FALSE,
  #15) data_type="Float32" # if NULL find it using dataType_table now
  #16) scaling=NULL,
  #17) values_range=NULL
  
  #data_type_str <- unique(dataType(test[[1]][[1]]$r_local_moran)
  r_test <- test[[1]][[1]]$r_local_moran
  plot(r_test)
  data_type_str <- unique(dataType(r_test))
                          
  NAvalue(r_test) #-inf hence use the dataType_table function
  range(r_test)
                          
  dataType_table <- generate_raster_dataType_table()
  
  dataType_selected <- dataType_table$r_type==data_type_str
  data_type_table_selected <- dataType_table[dataType_selected,]
  data_type_table_selected
  
  data_type <- data_type_table_selected$gdal_type
  
  min_val <- data_type_table_selected$min
  max_val <- data_type_table_selected$max
  valid_range <- c(min_val,max_val)
  NA_flag_val <- data_type_table_selected$min

  values_range <- NULL
  tmp_files <- FALSE #arg 18, param 18, keep temp files if TRUE
  scaling <- 1 #, param 20, if null use 1, for world mosaic use 1, since it is already multiplied by 100
  
  if(is.null(values_range)){
    values_range <- valid_range
  }

  mosaic_python <- "/nfs/bparmentier-data/Data/projects/spatial_variation_GIMMS/scripts"
  python_bin <- "/usr/bin/" #python gdal bin, on Atlas NCEAS
  match_extent <- "FALSE" #PARAM 31 #try without matching!!!
  #data_type <- NULL
  mosaicing_method <- "use_edge_weights" #PARAM10, arg 10
  algorithm <- "python" #PARAM 16 #if R use mosaic function for R, if python use modified gdalmerge script from Alberto Guzmann
  r_mask_raster_name <- NULL
  
  # list output tiles directories
  in_dir_tiles_tmp <- file.path(out_dir,paste0("tile_",1:n_tiles))
  
  pattern_str <- "lag_1_10.*.tif"
  ## get list of files for the mosaicing:
  lf <- lapply(in_dir_tiles_tmp,
               function(x){list.files(x,pattern=pattern_str,full.names = T)})
  n_dates <- length(lf[[1]])
  n_tiles <- length(in_dir_tiles_tmp)
  ## reorganize files for mosaicing:
  lf_mosaic <- get_files_to_mosaic(lf,n_dates,n_tiles)
  
  ### Use mean average?
  
  mosaic_method <- "use_edge_weights" #this is distance from edge
  day_to_mosaic <- paste("date_",1:n_dates,sep="")

  i <- 1
  #debug(mosaicFiles)
  
  list_date_mosaic_obj <- vector("list",length=length(lf_mosaic))
  
  ## For every date
  for(i in 1:length(lf_mosaic)){
    
    lf_in <- lf_mosaic[[i]]
    out_suffix_date <- paste(day_to_mosaic[i],mosaicing_method,out_suffix,sep="_")
    #out_dir_tmp <- paste0("date_outputs_date",i,"_tmp")
    out_dir_tmp <- paste0("date_outputs_",out_suffix_date)
    out_dir_tmp <- file.path(out_dir,out_dir_tmp)
    
    if(!dir.exists(out_dir_tmp)){
      dir.create(out_dir_tmp)
    }
    
    #debug(split_multiband)
    #lf_multi <- split_multiband(lf_in[1],out_suffix,out_dir_tmp)
    
    out_suffix_tmp <- "_temporary"
    lf_multi <- lapply(lf_in,
                       FUN=split_multiband,
                       out_suffix = out_suffix_tmp,
                       out_dir = out_dir_tmp)
      
    df_multiband <- as.data.frame(do.call(cbind,lf_multi))
    names(df_multiband) <- paste0("tile_",1:n_tiles)
    
    ### bug error, columns have become factor: changed this here
    df_multiband <- data.frame(lapply(df_multiband, as.character), stringsAsFactors=FALSE)

    j <- 1         
    ## for every band:
    list_mosaic_obj <- vector("list",length=nrow(df_multiband))
    
    for(j in 1:nrow(df_multiband)){
      ## problem with level variable: fixing this
      in_file <- as.list(df_multiband[j,])
      out_suffix_tmp_band <-  paste0("b_",j,"_",out_suffix_date)
      #out_suffix_date
      
      #undebug(mosaicFiles)
      
      ### All parameters should be set up earlier!!!!
      
      list_mosaic_obj[[j]] <- mosaicFiles(in_file,
                                mosaic_method="use_edge_weights",
                                num_cores=num_cores,
                                r_mask_raster_name=r_mask_raster_name, # can bin NULL
                                python_bin=python_bin,
                                mosaic_python=mosaic_python,
                                algorithm=algorithm,
                                match_extent=match_extent,
                                df_points=NULL,
                                NA_flag=NA_flag_val,
                                file_format=file_format,
                                out_suffix=out_suffix_tmp_band,
                                out_dir=out_dir_tmp,
                                tmp_files=tmp_files,
                                data_type=data_type,
                                scaling=scaling,
                                values_range=values_range)
    }
      
    #runs in 15-16 minutes for 3 dates and mosaicing of 28 tiles...
    list_date_mosaic_obj[[i]] <- list_mosaic_obj[[j]]
    
    #need to clean up temporary files here:
    if(tmp_files==F){ #if false...delete all files with "_tmp"
      lf_tmp <- list.files(pattern="*.*temporary*.*",path=out_dir_tmp,full.names=T)
      ##now delete temporary files...
      file.remove(lf_tmp)
    }
  }
    
  
}

#########################################
########### PART 5: Results: Examining lag information and variability ###########

#list.files(out_dir)
list_dirs_outputs <- list.dirs(out_dir)

list_dir_mosaics <- list_dirs_outputs[grepl("date_outputs_date.*.",list_dirs_outputs)]

lf <- mixedsort(list.files(list_dir_mosaics[[1]],
                           pattern="r_m_.*.mean_masked.*.tif",full.names=T))
r_stack <- stack(lf) 
#plot(r_stack)
plot(r_stack,y=2)

#debug(grid_sampling_raster)
grid_samples_sf <- grid_sampling_raster(x_sampling=4,y_sampling=4,r=r_stack)

grid_samples_sp <- as(grid_samples_sf,"Spatial")

x_val <- extract(r_stack,grid_samples_sp,df=T)
#class(x_val)
dim(x_val)
x_val <- as(x_val,"sf")
plot(as.numeric(x_val[1,2:11]),type="l")
plot(as.numeric(x_val[2,2:11]),type="l")
plot(as.numeric(x_val[3,2:11]),type="l")
plot(as.numeric(x_val[10,2:11]),type="l")
View(x_val)


animate(r_local_moran_stack) ## Generate movie later on:


#plot(x_val[1,],type="l",ylim=c(-1.2,1.2))
#lines(x_val[2,],type="l",col="red")
#lines(x_val[3,],type="l",col="green")
#lines(x_val[4,],type="l",col="blue")

names(moran_df) <- c("moran_I","lag")
plot(moran_df$moran_I ~moran_df$lag)
#plot(moran_df$moran_I ~moran_df$lag,ylim=c(0,1))

#prepare to automate the plotting   all columns
#mydata<-moran_df
#dd <- do.call(make.groups, mydata[,-ncol(mydata)]) 
#dd$lag <- mydata$lag 

#layout_m<-c(1,3) #one row three columns
#res_pix <-  500
#png(paste("Figure_0b_graphic_abstract_spatial_correlogram_tmax_prediction_models_gam_levelplot_",date_selected,out_prefix,".png", sep=""),
#    height=res_pix*layout_m[1],width=res_pix*layout_m[2])

#p<-xyplot(data ~ lag | which, dd,type="b",main="Spatial content in interpolated Surfaces using GAM on September 1, 2010",
#          par.settings = list(axis.text = list(font = 2, cex = 1.3),layout=layout_m,
#                              par.main.text=list(font=2,cex=2),strip.background=list(col="white")),par.strip.text=list(font=2,cex=1.5),
#          strip=strip.custom(factor.levels=names_layers),
#          xlab=list(label="Spatial lag neighbor", cex=2,font=2),
#          ylab=list(label="Moran's I", cex=2, font=2))
#print(p)

#dev.off()

##########################################  END OF SCRIPT  ##############################
