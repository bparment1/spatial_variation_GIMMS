############### SESYNC Research Support: GIMMS Spatial Variation ########## 
## Help in processing modis data for the workshop group.
##
##
## DATE CREATED: 01/24/2018
## DATE MODIFIED: 06/27/2018
## AUTHORS: Benoit Parmentier  
## Version: 1
## PROJECT: spatial variability landscape
## ISSUE: 
## TO DO:
##
## COMMIT: Splitting functions and main scripts
##
## Links to investigate:
#https://cran.r-project.org/web/packages/gstat/vignettes/st.pdf

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

raster_processing_functions <- "GIMMS_processing_functions_02062018.R" #Functions used to mosaic predicted tiles
generate_tiles_functions <- "generate_spatial_tiles_functions_06282018.R" #Functions used to mosaic predicted tiles
source(file.path(script_path,raster_processing_functions)) #source all functions used in this script 
source(file.path(script_path,generate_tiles_functions))

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
out_suffix <-"GIMMS_processing_06272018" #output suffix for the files and ouptut folder
#ARGS 8
num_cores <- 2 # number of cores
#ARGS 9
date_param <- "1982.01.01;1982.12.31" #start date, end date
#ARGS 10
GIMMS_product <- "3g.v1"
#ARGS 11
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
  f <- list.files(path=in_dir,pattern="*.nc4",full.names=T)
}


##############################
#### PART 2: Import dataset
####
##### Next import in Tif format the NCDF

if(processing_steps$download==TRUE){
  
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
  
}

#### Next steps to consider:
## Use QC index to screen for low value pixels

#########################################
########### PART 3: Generate grid ###########

#This can be non-overlapping or overlapping

lf_gimms <- mixedsort(list.files(pattern=file_format,path=in_dir,full.names=T))

ref_file <- lf_gimms[1]

##### Generate a grid/tile for processing:
## Must transformed to a function later on.

if(processing_steps$tiling==TRUE){
  
  r_ref <- raster(ref_file)
  
  plot(r_ref)
  
  ### generate overlapping tiles
  x_overlap <- 0.25
  y_overlap <- 0.25
  
  ### store in the tile directory
  #undebug(generate_tiles_from_extent)
  
  #we want 20x20
  x_ratio <- 360/20
  y_ratio <- 180/20
  
  test_tiles2_2_overlap50 <- generate_tiles_from_extent(r_ref,
                                                        y_ratio=y_ratio,
                                                        x_ratio=x_ratio,
                                                        y_overlap=y_overlap,
                                                        x_overlap=x_overlap,
                                                        out_suffix=out_suffix,
                                                        out_dir=NULL) #this does not work...
  
  df_tiles <- test_tiles2_2_overlap50$df_tiles
  
  
  
  head(df_tiles)
  dim(df_tiles) #too many tiles: 288
  
  plot(r_ref)
  plot(test_tiles2_2_overlap50$list_tiles[[89]],add=T)
  
  #### Transform a data.frame object to a sf object using coordinates x and y
  centroids_tiles <- st_as_sf(df_tiles,
                      coords = c('x_center', 'y_center'),
                      crs = proj4string(r_ref))
  
  #plot(centroids_tiles["ID"],add=T)
  text(df_tiles$x_center,df_tiles$y_center,df_tiles$ID,cex=0.5)
  centroids_tiles$ID
  
  out_centroids_filename <- file.path(out_dir,"tiles_centroids.shp")
  st_write(centroids_tiles,
           dsn=out_centroids_filename)
  
  #Select tiles matching the areas of interest and with specific content:
  #Kenya<-getData("GADM", country="KE", level=0)
  KZ_spdf <-getData("GADM", country="KZ", level=0)
  KZ_sf <- as(KZ_spdf,"sf")
  
  plot(KZ_sf$geometry,add=T,border="blue")
  #Kenya1<-getData("GADM", country="KE", level=1)
  
  tiles_combined_spdf <- do.call(rbind,test_tiles2_2_overlap50$list_tiles)
  plot(tiles_combined_spdf)
  tiles_combined_sf <- as(tiles_combined_spdf,"sf") #transforming spdf to sf object
  st_crs(KZ_sf)
  st_crs(tiles_combined_sf) <- proj4string(r_ref)
  
  tiles_combined_sf$ID <- 1:nrow(tiles_combined_sf)
  tiles_all <- join(tiles_combined_sf,df_tiles)
  tiles_combined_sf <- merge(tiles_combined_sf,df_tiles,by="ID")
  View(tiles_combined_sf)
  
  test_sf <- st_intersection(tiles_combined_sf,KZ_sf)
  test_sf <- st_intersection(KZ_sf,tiles_combined_sf)
  test_sf
  plot(test_sf$geometry,border="red")
  #test_sf <- st_intersects(tiles_combined_sf,KZ_sf)
  test_sf <- st_intersects(KZ_sf,tiles_combined_sf)
  test_sf[[1]]
  class(test_sf)
  unlist(test_sf)
  tiles_combined_sf$ID
  selected_tiles <- filter(tiles_combined_sf,ID%in%unlist(test_sf))
  #selected_tiles
  
  #sum(tiles_combined_sf$ID==unlist(test_sf))
  
  
  plot(r_ref)
  text(df_tiles$x_center,df_tiles$y_center,df_tiles$ID,cex=0.5)
  plot(KZ_sf$geometry,add=T,border="blue")
  #plot(test_tiles2_2_overlap50$list_tiles[[89]],add=T)
  plot(selected_tiles$geometry,border="red",add=T)
  
  ## Let's drop 41 isnce the overlap is low
  
  selected_tiles
  
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

#grid_filename
grid_filename <- out_grid_filename
tile_grid <- st_read(grid_filename)

tile_index <- 56

### Plot tiles being processed:

#tile_grid_selected <- as(tile_grid[tile_index,],"Spatial")
tile_grid_selected <- tile_grid[tile_index,]

tile_selected_spdf <- as(tile_grid_selected, "Spatial")

plot(r)
plot(tile_grid$geometry,add=T)
plot(tile_grid_selected,add=T,col=NA,border="red")

text(gCentroid(tile_selected_spdf),paste0("tile ",tile_index, " processed"))
#text(st_centroid(tile_grid_selected),paste0("tile ",tile_index, " processed"))
r_tile <- crop(r,tile_selected_spdf)

#writeRaster(r_tile,)

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

list_filters<-lapply(1:10,FUN=autocor_filter_fun,f_type="queen") #generate 10 filters
#moran_list <- lapply(list_filters,FUN=Moran,x=r)

r_stack <- r_tile

list_param_moran <- list(list_filters=list_filters,
                         r_stack=r_stack,
                         out_suffix=NULL,
                         out_dir=NULL)

#moran_r <-moran_multiple_fun(1,list_param=list_param_moran)
nlayers(r_stack) 

strsplit(names(r),split="[.]")

#debug(local_moran_multiple_fun)
r_test <- local_moran_multiple_fun(1,list_param=list_param_moran)

  
local_moran_I_list <-mclapply(1:nlayers(r_stack), list_param=list_param_moran, 
                      FUN=local_moran_multiple_fun,mc.preschedule=FALSE,
                      mc.cores = 2) #This is the end bracket from mclapply(...) statement

r_local_moran_stack <- stack(unlist(local_moran_I_list))

#########################################
########### PART 5: Results: Examining lag information and variability ###########

plot(r_local_moran_stack,y=2)
animate(r_local_moran_stack)

#### Generate longitude and latitude info
locations_mat <- rbind(c(-110.6982,34.4208),
                      c(-119.6982,34.4208),
                      c(-101,42),
                      c(-123,47))

locations_df <- as.data.frame(locations_mat)
names(locations_df)<- c("longitude","latitude")
locations_df$ID <- 1:nrow(locations_df)
locations_sf <- st_as_sf(locations_df,coords=c("longitude", "latitude"), 
                         crs = 4326) 
  
  
plot(r_tile)
plot(locations_sf$geometry,add=T)

x_val <- extract(r_local_moran_stack,locations_mat)

plot(x_val[1,],type="l",ylim=c(-1.2,1.2))
lines(x_val[2,],type="l",col="red")
lines(x_val[3,],type="l",col="green")
lines(x_val[4,],type="l",col="blue")


moran_I_df <-mclapply(1:nlayers(r_stack), list_param=list_param_moran, 
                      FUN=moran_multiple_fun,mc.preschedule=FALSE,
                     mc.cores = 2) #This is the end bracket from mclapply(...) statement

moran_df <- do.call(cbind,moran_I_df) #bind Moran's I value 10*nlayers data.frame
moran_df$lag <-1:nrow(moran_df)

names(moran_df) <- c("moran_I","lag")
plot(moran_df$moran_I ~moran_df$lag)
#plot(moran_df$moran_I ~moran_df$lag,ylim=c(0,1))

#prepare to automate the plotting of   all columns
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
