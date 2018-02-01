############### SESYNC Research Support: Hurricane Management ########## 
## Help in processing modis data for the workshop group.
##
##
## DATE CREATED: 01/24/2018
## DATE MODIFIED: 02/02/2018
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

###### Functions used in this script

#Benoit setup
script_path <- "/nfs/bparmentier-data/Data/projects/spatial_variation_GIMMS/scripts"

raster_processing_functions <- "GIMMS_processing_functions_02012018.R" #Functions used to mosaic predicted tiles
source(file.path(script_path,raster_processing_functions)) #source all functions used in this script 

############################################################################
#####  Parameters and argument set up ###########

#ARGS 1
#Set up teamhurricane-data
in_dir <- "/nfs/bparmentier-data/Data/projects/spatial_variation_GIMMS/data"

out_dir <- "/nfs/bparmentier-data/Data/projects/spatial_variation_GIMMS/outputs"

#hdf_file <-"test.hdf"

#NA_flag <- -999999
NA_flag_val <- NULL

file_format <- ".tif"
scaling_factor <- 0.0001 #MODIFY THE SCALING FACTOR - FOR NORMALIZED DATA SHOULD BE 10,000 AT LEAST
#ARGS 7
create_out_dir_param=TRUE #create a new ouput dir if TRUE
#ARGS 8
out_suffix <-"GIMMS_processing_02012018" #output suffix for the files and ouptut folder #param 12
num_cores <- 2 # number of cores

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
### PART I DOWNLOAD AND PREPARE DATA #######

#Will be a function

## Information: https://nex.nasa.gov/nex/projects/1349/wiki/general_data_description_and_access/

## Download data here:
#https://ecocast.arc.nasa.gov/data/pub/gimms/3g.v1/
#lf <- list.files(dir="https://ecocast.arc.nasa.gov/data/pub/gimms/3g.v1/")

lf_name <- download.file("https://ecocast.arc.nasa.gov/data/pub/gimms/3g.v1/00FILE-LIST.txt",
                         destfile = "00FILE-LIST.txt")

lf_df <- read.table("00FILE-LIST.txt",stringsAsFactors = F)

## Make this a function later on!!!

date_param <- "2002.01.01;2012.12.31;8" #start date, end date, time_step

start_date <- date_param[1]
end_date <- date_param[2]

GIMMS_product <- "3g.v1"

debug(modis_product_download)
modis_product_download(GIMMS_product,
                       start_date,
                       end_date,
                       out_dir,
                       out_suffix)  ##Functions used in the script
  

##### Next import in Tif format the NCDF

GDALinfo_raster <- GDALinfo(raster_file,returnScaleOffset = F) #use GDAL info utility

#hdf_file <- file.path(in_dir,hdf_file) #associate file path to input
#GDALinfo_hdf <- GDALinfo(hdf_file,returnScaleOffset = F) #use GDAL info utility

#str(GDALinfo_hdf) ## Class GDLobj
raster_subdataset <- attributes(GDALinfo_raster)$subdsmdata
print(raster_subdataset)

### Generate a data.frame from scientific datasets
raster_df <- (strsplit(raster_subdataset,":"))
length(raster_subdataset)
raster_df <- as.data.frame(do.call(rbind,raster_df),stringsAsFactors=F)
raster_df <- data.frame(lapply(raster_df, as.character))

#hdf_df %>% 
#  mutate_all(as.character)

#names(raster_df) <- c("subdataset_name","description","dir","product","var_name")
names(raster_df) <- c("subdataset_name","dir","var_name")

raster_df

#Select automatically QC flag!!
#View(hdf_df)

write.table(raster_df,"raster_subdataset.txt",sep=",")

#modis_subset_layer <- paste("HDF4_EOS:EOS_GRID:",
#                                raster_file,subdataset,sep="")

#raster_subset_layer <- paste0("NETCDF:","ndvi3g_geo_v1_1981_0712.nc4",":","ndvi")
#r <- readGDAL(raster_subset_layer) #read specific dataset in hdf file and make SpatialGridDataFrame
#r  <- brick(r) #convert to raser object

lf<- mixedsort(list.files(pattern="*.nc4"))

#lf[2]

#### New import function can be created here...
i <- 2

raster_file <- lf[i]
r <- brick(raster_file,"ndvi")

raster_name <- sub(extension(raster_file),"",raster_file)
list_raster_name <- unlist(strsplit(x=basename(raster_name), split="[_]"))
year_val <- list_raster_name[[length(list_raster_name)-1]]
month_range <- list_raster_name[[length(list_raster_name)]]
month_range <- c(substr(month_range, 1, 2),substr(month_range, 3, 4))
#month_range <- as.numeric(month_range)

start_date <- paste(year_val,month_range[1],"01",sep=".")
end_date <- paste(year_val,month_range[2],"30",sep=".")

st <- as.Date(start_date,format="%Y.%m.%d") #start date
en <- as.Date(end_date,format="%Y.%m.%d") #end date
ll <- seq.Date(st, en, by="15 day") #sequence of dates
#dates_queried <- format(ll,"%Y.%m.%d") #formatting queried dates
dates_val <- format(ll[-13],"%Y_%m_%d") #formatting queried dates


#generate dates for 16 days product
#dates_val <- generate_dates_by_step(date_range[1],date_range[2],16)$dates #NDVI Katrina
#names_ncdf <- as.character(unlist(strsplit(x=basename(raster_file), split="[_]")))
#names_ncdf <- names_ncdf[-length(names_ncdf)]
#names_ncdf <- names_ncdf[-length(names_ncdf)]
#strsplit(names_ncdf,split="[.]")

##Assign new names for r?
#names_layers <- names(r)
#names_layers <- gsub("X","m_",names_layers)

#names(r)
names(r) <- paste("NDVI",dates_val,sep="_")
#now drop "." and also later right out the date in the file!!!

raster_name <- paste0(paste(names_ncdf,collapse = "_"),"_",out_suffix,file_format)


#gsub("/.","_",names_layers)
#sub(".","_",names_layers)
#lapply(names_layers,function(x){sub(".","_",x)})
#names(r) <-

if(is.null(NA_flag_val)){
  
  NA_flag_val <- NAvalue(r)
  
}

#NDVI variable
#modis_layer_str1 <- unlist(strsplit(modis_subdataset[1],"\""))[3] #Get day NDVI layer
#QC
#modis_layer_str2 <- unlist(strsplit(modis_subdataset[5],"\""))[3] #Get day VI QC layer

#subdataset <- modis_layer_str1
#modis_subset_layer_Day <- paste("HDF4_EOS:EOS_GRID:",hdf_file,subdataset,sep="")

#r <- readGDAL(modis_subset_layer_Day) #read specific dataset in hdf file and make SpatialGridDataFrame

plot(r,y=1,zlim=c(-10000,10000))
r #print propind the NA flag value
data_type_str <- dataType(r) #find the dataTyre

# for more control, you can set dataType and/or compress the files
#data_type_str <- "FLT4S"
#A_flag_val <- NAvalue(r2)

#This is a brick!!!! so need to write out every layer!!

#raster_name <- ""
  
#rwriteRaster(r,
#            file.path(out_dir,raster_name),
#            overwrite=TRUE,
#            NAflag=NA_flag_val,
#            datatype=data_type_str,
#            options=c("COMPRESS=LZW"))

#writeRaster(r_stack,
#            filename=file.path(out_dir,paste("pysal_mle",".rst",sep="")),
#            bylayer=T,suffix=paste(names(r_stack),"_",out_suffix,sep=""),
#            overwrite=TRUE,
#            NAflag=NA_flag_val,
#            datatype=data_type_str,
#            options=c("COMPRESS=LZW"))

#names_hdf<-as.character(unlist(strsplit(x=basename(hdf), split="[.]")))

names_ncdf <-as.character(unlist(strsplit(x=basename(raster_file), split="[_]")))
names_ncdf <- names_ncdf[-length(names_ncdf)]
names_ncdf <- names_ncdf[-length(names_ncdf)]
strsplit(names_ncdf,split="[.]")


#raster_name <- 
out_suffix_str <- ""
#raster_name <- paste0(paste(names_ncdf,collapse = "_"),"_",out_suffix,file_format)
raster_name <- paste0(paste(names_ncdf,collapse = "_"),file_format)

if(file_format==".tif"){
  format_raster="GTiff"
}

writeRaster(r,
            filename=file.path(out_dir,raster_name,sep=""),
            bylayer=T,
            #suffix=paste(names(r),"_",out_suffix,sep=""),
            format=format_raster,
            suffix=paste(names(r)),
            overwrite=TRUE,
            NAflag=NA_flag_val,
            datatype=data_type_str,
            options=c("COMPRESS=LZW"))

#### Next steps to consider:
## Use the name from MODIS file because it contains information on tile location, date and product type
## Use QC index to screen for low value pixels

########### PART 3: Now do the analyses ###########

#Tthis can be in another script.

lf_gimms <- mixedsort(list.files(pattern=file_format,path="."))

r <- raster(lf_gimms[1])

##### Generate a grid/tile for processing:
## Must transformed to a function later on.

# for the time being generate a no-overlapping grid tiling and crop
extent_val <- extent(r)
bbox_val <- st_bbox(r)
test_sp <- as(extent_val, 'SpatialPolygons')
outline_sf <-as(test_sp,"sf")

#Can buffer?

#test_grid <- st_make_grid(outline_sf, n=18)
test_grid <- st_make_grid(outline_sf, n=9)

plot(r)
plot(test_grid,add=T)
plot(test_grid[56],add=T,col="red")
tile_grid_selected <- as(test_grid[56],"Spatial")
r_tile <- crop(r,tile_grid_selected)
  
### Create function for SLURM job:


#ARGS 11
infile_list_tiles <- "list_tiles.txt"

#ARGS 12: 
#TILE INDEX this is from the array
# Get the the array id value from the environment variable passed from sbatch

SLURM_ARRAY_TASK_ID <- Sys.getenv('SLURM_ARRAY_TASK_ID')
tile_index <- Sys.getenv("SLURM_ARRAY_TASK_ID") #this is should be an integer from 1:n, n is the number of tiles
#slurm_arrayid <- Sys.getenv('SLURM_ARRAYID') #work with #SARRAY option in SLURM
#tile_index <- as.numeric(slurm_arrayid) # coerce the value to an integer
#tile_index <- 1  #for testing

#generate filters for 10 lags: quick solution

list_filters<-lapply(1:10,FUN=autocor_filter_fun,f_type="queen") #generate 10 filters
#moran_list <- lapply(list_filters,FUN=Moran,x=r)

r_stack <- r_tile

##
list_param_moran <- list(list_filters=list_filters,
                         r_stack=r_stack,
                         out_suffix=NULL,
                         out_dir=NULL)

#moran_r <-moran_multiple_fun(1,list_param=list_param_moran)
nlayers(r_stack) 

strsplit(names(r),split="[.]")

debug(local_moran_multiple_fun)
r_test <- local_moran_multiple_fun(1,list_param=list_param_moran)

  
local_moran_I_list <-mclapply(1:nlayers(r_stack), list_param=list_param_moran, 
                      FUN=local_moran_multiple_fun,mc.preschedule=FALSE,
                      mc.cores = 2) #This is the end bracket from mclapply(...) statement

r_local_moran_stack <- stack(unlist(local_moran_I_list))
plot(r_local_moran_stack,y=2)
animate(r_local_moran_stack)

x_val <- extract(r_local_moran_stack,cbind(-119.6982,34.4208))

plot(x_val[1,],type="l")

x_val <- extract(r_local_moran_stack,cbind(-110.6982,34.4208))
plot(x_val[1,],type="l")

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

######################### END OF SCRIPT ##############################
