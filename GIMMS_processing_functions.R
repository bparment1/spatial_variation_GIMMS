############### SESYNC Research Support: Hurricane Management ########## 
## Help in processing modis data for the workshop group.
##
##
## DATE CREATED: 01/29/2018
## DATE MODIFIED: 02/02/2018
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

extractFolders=function(urlString) {
  htmlString=getURL(urlString)
  ret=gsub("]", "", str_replace_all(str_extract_all(htmlString, paste('DIR',".([^]]+).", '/\">',sep=""))[[1]], "[a-zA-Z\"= <>/]", ""))
  return(ret[which(nchar(ret)>0)])
}

#list_folders_files[[i]] <- extractFiles(url_folders_str[i], list_tiles)[file_format]
extractFiles=function(urlString, list_tiles_str) {
  #Slight modifications by Benoit
  #list_tiles: modis tiles as character vectors eg c("h10v06","h09v07")
  #urlString: character vector with url folder to specific dates for product
  
  # get filename strings
  htmlString=getURL(urlString)
  #htmlString=getURL(urlString[2])
  allVec=gsub('\">', '', gsub('<a href=\"', "", str_extract_all(htmlString, paste('<a href=\"',"([^]]+)", '\">',sep=""))[[1]]))
  #allVec: this contains list of all files! need to select the correct tiles...
  #ret=c()
  #for (currSel in list_tiles_str) {
  #  ret=c(ret, grep(currSel, allVec, value=TRUE))
  #}
  #list_tiles_str <- c("h10v06","h11v07")
  ret <- lapply(list_tiles_str,function(x){grep(x,allVec,value=T)})
  
  # select specific files
  #ret <- paste(urlString,ret,sep="") #append the url of folder
  ret <-file.path(urlString,unlist(ret))
  jpg=sapply(ret, FUN=endswith, char=".jpg")
  xml=sapply(ret, FUN=endswith, char=".xml")
  hdf=sapply(ret, FUN=endswith, char=".hdf")
  
  retList=list(jpg=ret[which(jpg)], xml=ret[which(xml)], hdf=ret[which(hdf)])
  return(retList)
}

endswith=function(x, char) {
  currSub = substr(x, as.numeric(nchar(x)-nchar(char))+1,nchar(x))
  if (currSub==char) {return(TRUE)}
  return(FALSE)
}

GIMMS_product_download <- function(GIMMS_product,start_date,end_date,out_dir,out_suffix){  
  
  ########## BEGIN SCRIPT #######
  
  #step 1: parse input elements
  
  st <- as.Date(start_date,format="%Y.%m.%d") #start date
  en <- as.Date(end_date,format="%Y.%m.%d") #end date
  
  year(st)
  year(en)
  
  #date_param <- "2002.01.01;2012.12.31;8" #start date, end date, time_step
  
  lf_name <- download.file("https://ecocast.arc.nasa.gov/data/pub/gimms/3g.v1/00FILE-LIST.txt",
                           destfile = "00FILE-LIST.txt")
  
  lf_df <- read.table("00FILE-LIST.txt",stringsAsFactors = F)
  
  ## Make this a function later on!!!
  
  nf <- 3 #number of files to download
  list_raster_file <- vector("list",length=nf)
  for(i in 1:nf){
    raster_file <- basename(lf_df[i,1]) #this is the outfile
    
    file1 <- download.file(lf_df[i,1],raster_file)
    #https://ecocast.arc.nasa.gov/data/pub/gimms/
    
    list_raster_file[[i]] <- raster_file
  }
  

  #Prepare return object: list of files downloaded with http and list downloaded of files in tiles directories
  
  list_files_by_tiles <-mapply(1:length(out_dir_tiles),
                               FUN=function(i,x){list.files(path=x[[i]],pattern="*.hdf$",full.names=T)},MoreArgs=(list(x=out_dir_tiles))) #Use mapply to pass multiple arguments
  #list_files_by_tiles <-mapply(1:length(out_dir_tiles),FUN=list.files,MoreArgs=list(pattern="*.hdf$",path=out_dir_tiles,full.names=T)) #Use mapply to pass multiple arguments
  
  colnames(list_files_by_tiles) <- list_tiles #note that the output of mapply is a matrix
  download_modis_obj <- list(list_files_tiles,list_files_by_tiles)
  names(download_modis_obj) <- c("downloaded_files","list_files_by_tiles")
 
   return(download_modis_obj)
}


import_gimms_nc4 <- function(f,out_dir,var_name="ndvi",out_suffix="",out_dir="."){
  
  #var_name <- "ndvi"
  
  #### Begin script ###
  
  #ll <- seq.Date(st, en, by="1 day") #sequence of dates
  #dates_queried <- format(ll,"%Y.%m.%d") #formatting queried dates
  
  #year(date_queried)
  #raster_file <- lf[i]
  raster_file <- f
  
  r <- brick(raster_file,var_anem)
  
  raster_name <- sub(extension(raster_file),"",raster_file)
  list_raster_name <- unlist(strsplit(x=basename(raster_name), split="[_]"))
  year_val <- list_raster_name[[length(list_raster_name)-1]]
  month_range <- list_raster_name[[length(list_raster_name)]]
  month_range <- c(substr(month_range, 1, 2),substr(month_range, 3, 4))
  #month_range <- as.numeric(month_range)
  
  #start_date <- paste(year_val,month_range[1],"01",sep=".")
  #end_date <- paste(year_val,month_range[2],"30",sep=".")
  
  #st <- as.Date(start_date,format="%Y.%m.%d") #start date
  #en <- as.Date(end_date,format="%Y.%m.%d") #end date
  #ll <- seq.Date(st, en, by="15 day") #sequence of dates
  
  month_vals <- month_range[1]:month_range[2]
  month_vals <- sprintf("%02d", month_vals)
  
  list_dates <-unlist(lapply(month_vals,function(x){paste0(x,".",c("01","15"))}))
  #dates_queried <- format(ll,"%Y.%m.%d") #formatting queried dates
  #dates_val <- format(ll[-13],"%Y_%m_%d") #formatting queried dates
  
  
  
  
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
  
}

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
  #1) list_filters: list of filters with different lags in the image
  #2) r_stack: stack of raster image, only the selected layer is used...
  #3) out_suffix: if NULL default, no suffix added.
  #4) out_dir: if NULL default, use current dir
  #OUTPUTS:
  #
  
  ###### Begin #######
  
  if(is.null(out_suffix)){
    out_suffix <- ""
  }
  if(is.null(out_dir)){
    out_dir <- "."
  }
  
  ### Read in filters used to computer Moran's I
  
  list_filters <-list_param$list_filters
  
  r <- subset(list_param$r_stack,i)
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
  
  writeRaster(r_local_moran,
              filename=raster_name,
              bylayer=T,
              suffix=paste(names(r_local_moran),"_",out_suffix,sep=""),
              overwrite=TRUE,
              NAflag=NA_flag_val,
              datatype=data_type_str,
              options=c("COMPRESS=LZW"))
  
  #### Prepare return object
  
  return(moran_list)
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

