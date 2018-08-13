##############################  Raster Processing Workflow  #######################################
#######################  Spatial utility: General code for mosaicing tiles/spatial subsets  ##############################
#This script contains utility functions to mosaics tiles for given areas or the globe.
#Different options to explore mosaicing are tested. This script only contains functions.
#AUTHOR: Benoit Parmentier 
#CREATED ON: 04/14/2015  
#MODIFIED ON: 08/14/2018            
#Version: 2
#PROJECT: Environmental Layers project     
#COMMENTS:
#
#COMMIT: deal with raster temporary files from raster package in mosaics
#
#TODO:
#1) Make this is a script/function callable from the shell/bash
#2) Improve performance: there will be a need to improve efficiency for the workflow.

#Error message for gdal_proximity:
#ERROR 1: Source and proximity bands are not the same size.
#gdal_proximity.py give an error when tries to replace an existent output file with different band.
#If the output already exists and was created by a different band input file, this error is displayed: ERROR 1: Source and proximity bands are not the same size.
#If you remove the existent file, it works fine.
#available:
#See below

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

#### FUNCTION USED IN SCRIPT

##List all the functions in this script:

#[2] "create_accuracy_metric_raster"     
#[3] "create_accuracy_residuals_raster"   
#[4] "create_weights_fun"                
# [5] "fit_models"                         "function_mosaicing"                
# [7] "in_dir_script"                      "mosaicFiles"                       
# [9] "mosaic_m_raster_list"               "mosaic_python_merge"               
#[11] "plot_daily_mosaics"                 "plot_diff_raster"                  
#[13] "plot_mosaic"                        "plot_screen_raster_val"            
#[15] "predict_accuracy_raster_by_station" "predict_auto_krige_raster_model"   
#[17] "raster_match"                       "remove_na_spdf"                    
#[19] "select_var_stack"                   "sine_structure_fun"    

############## START FUNCTIONS DEFINITIONS ####

sine_structure_fun <-function(x,T,phase,a,b,use_cos=FALSE){
  
  #Create sine for a one dimensional series
  #Note that sine function uses radian unit.
  #a=amplitude
  #b=mean or amplitude 0 of the series
  #T= stands for period definition
  #phase=phase angle (in radian!!)
  #cos: use cosine instead of sine if TRUE
  
  if(use_cos==FALSE){
    y <- a*sin((x*pi/T)+ phase) + b
  }else{
    y <- a*cos((x*pi/T)+ phase) + b
  }
  return(y)
}

create_raster_df_centroids_fun <- function(j,list_param){
  #This function generates raster images from metrics values from a data.frame.
  #The raster layer is assigned a unique value from the pixel at every location.
  #Input Parameters:
  #  #lf: list of raster files
  #df_centroids: data.frame table #fitting or validation table with all days
  #metric_name: accuracy metric selected to be mapped, RMSE, MAE etc.
  #date_processed: day being processed , e.g. 19920101
  #num_cores : number of cores used in the parallelization
  #NA_flag_val: value used as flag in the raster 
  #file_format: e.g. tif., .rst
  #out_dir_str: output directory
  #out_suffix_str: output suffix
  #Outputs:
  
  #### PARSE arguments
  
  df_centroids <- list_param$df_centroids
  metric_name <- list_param$metric_name 
  #interpolation_method <- list_param$interpolation_method #c("gam_CAI") #PARAM3
  #date_processed <- list_param$metric_name 
  #num_cores <- list_param$num_cores
  NA_flag_val <- list_param$NA_flag_val
  file_format <- list_param$file_format
  out_dir_str <- list_param$out_dir
  out_suffix_str <- list_param$out_suffix
  
  ####### START SCRIPT #####    
  
  inFilename <- df_centroids$files[j]
  r1 <- raster(inFilename)
  r1[] <- df_centroids[[metric_name]][j] #improve this
  #set1f <- function(x){rep(NA, x)}
  #r_init <- init(r_in, fun=set1f)
  lf_tmp <- gsub(file_format,"",lf)
  
  extension_str <- extension(inFilename)
  raster_name_tmp <- gsub(extension_str,"",basename(inFilename))
  outFilename <- file.path(out_dir_str,paste(raster_name_tmp,"_",metric_name,"_",out_suffix,file_format,sep="")) #for use in function later...
  
  writeRaster(r1, NAflag=NA_flag_val,filename=outFilename,overwrite=TRUE)  
  #list_raster_name[[j]] <- outFilename
  return(outFilename)
}


mosaic_python_merge <- function(NA_flag_val,module_path,module_name,input_file,out_mosaic_name,raster_ref_name=NULL){
  #Inputs:
  #NA_falg_val: NA value to use in the output
  #module_path: python module to be used with path
  #module_name: name of the module
  #input_file: input file as text file containing a list of files to mosaic
  #out_mosaic_name: output file for the mosaic
  #rast_ref_name: reference file with extent and resolution to be used in the mosaicing
  #
  
  #### Start script ###
  
  
  if(is.null(raster_ref_name)){
    
    #out_mosaic_name <- r_weights_sum_raster_name <- file.path(out_dir,paste("r_weights_sum_m_",method_str,"_weighted_mean_",out_suffix,".tif",sep=""))
    cmd_str <- paste("python", file.path(module_path,module_name),
                     "--config GDAL_CACHEMAX=1500",
                     "--overwrite=TRUE",
                     paste("-o",out_mosaic_name,sep=" "),
                     paste("--optfile", input_file,sep=" "),
                     paste("-n",NA_flag_val,sep=" "))
    system(cmd_str)
  }
  if(!is.null(raster_ref_name)){
    
    #lf_files <- c(r_m_weighted_mean_raster_name) #match to mask
    #rast_ref_name <- r_mask_raster_name
    r_ref <- raster(raster_ref_name)
    extent_r_ref <- as.numeric(as.matrix(extent(r_ref)))
    res_pix <- res(r_ref)
    #c(xmin,ymax,xmax,ymin)
    #c(ulx,uly,lrx,lry)
    extent_str <- c(extent_r_ref[1],extent_r_ref[4],extent_r_ref[3],extent_r_ref[2])
    #-ps 0.00833349 0.008333229 -ul_lr -187.143383752 81.744912045 -5.073165638 14.554089652
    #out_mosaic_name <- r_weights_sum_raster_name <- file.path(out_dir,paste("r_weights_sum_m_",method_str,"_weighted_mean_",out_suffix,".tif",sep=""))
    cmd_str <- paste("python", file.path(module_path,module_name),
                     "--config GDAL_CACHEMAX=1500",
                     "--overwrite=TRUE",
                     paste("-o",out_mosaic_name,sep=" "),
                     paste("-ps", res_pix[1],res_pix[2],sep=" "), #pixel size
                     paste("-ul_lr", extent_str[1],extent_str[2],extent_str[3],extent_str[4],sep=" "), #extent
                     paste("--optfile", input_file,sep=" "),
                     paste("-n",NA_flag_val,sep=" "))
    system(cmd_str)
    
  }
  #list(out_mosaic_name,cmd_str)
  mosaic_python_merge_obj <- list(out_mosaic_name,cmd_str)
  names(mosaic_python_merge_obj) <- c("out_mosaic_name","cmd_str")
  
  return(mosaic_python_merge_obj)
}


create_weights_fun <- function(i, list_param){
  #This function generates weights from a point location on a raster layer.
  #Weights rasters are defined using a distance from edges or other option.s
  #Note that the weights are normatlized on 0-1 scale using max and min values.
  #Once the weights are rescaled, these are multiplied with the original 
  #input values from the raster to produce product weights by pixels.
  #There are three options for the weights:
  # - "use_sine_weights"
  # -  "use_edge"
  # - "use_linear_weights"
  #
  #INPUTS:
  #1)lf: list of raster files
  #2)df_points: reference points from which to compute distance
  #3)r_feature: reference features as raster image from which to compute distance from
  #4)methods: options available: use_sine_weights,use_edge,use_linear_weights
  #5)NA_flag : raster flag values, e.g. -9999
  #6)file_format: raster format used, default is ".tif"
  #7)out_suffix_str: output suffix, default is NULL, it is recommended to add the variable name etc.
  #             here e.g. dailyTmax and date!!
  #8)out_dir_str: output directory, default is NULL
  
  #OUTPUTS:
  #raster list of weights raster and product of weights raster and inputs
  #TODO: 
  # -use gdal proximity for large files and use_edge option
  # - add raster options
  # - improve efficiency
  # - change name options
  #
  ############
  
  ##### START SCRIPT #####
  
  ##### Parse out the input parameters
  
  lf <- list_param$lf
  df_points <- list_param$df_points
  r_feature <- list_param$r_feature #this should be change to a list
  padding <- TRUE #if padding true then make buffer around edges??
  method <- list_param$method #differnt methods available to create weights
  #NAflag,file_format,out_suffix etc...
  NA_flag_val <- list_param$NA_flag
  file_format <- list_param$file_format
  out_suffix_str <- list_param$out_suffix_str
  out_dir_str <- list_param$out_dir_str
  
  ##### Prepare weight layers  
  #lf <- unlist(lf)
  
  r_in <- raster(lf[i]) #input image
  tile_no <- i #file being processed, assuming tiles by tiles
  
  set1f <- function(x){rep(NA, x)}
  r_init <- init(r_in, fun=set1f)
  
  if(!is.null(r_feature)){
    r_init <- r_feature
  }
  
  if(!is.null(df_points)){ #reference points as SPDF object
    cell_ID <- cellFromXY(r_init,xy=df_points[i,])
    r_init[cell_ID] <- df_points$ID[i]
  }
  
  if(method=="use_sine_weights"){
    #Generate spatial pattern 5:     
    n_col <- ncol(r_init)
    n_row <- nrow(r_init)
    
    #u <- xFromCol(r_init,col=1:n_col)
    #add padding option later...buffer from a specific distance and tailling of at 0.1
    u <- 1:n_col
    a<- 1 #amplitude in this case
    b<- 0
    T<- n_col
    phase <- 0
    use_cos <- FALSE
    ux <- sine_structure_fun(u,T,phase,a,b,use_cos)
    ux_rep <-rep(ux,time=n_row)  
    r1 <-setValues(r_init,ux_rep)  #note efficient in memory might need to revise this
    #plot(r)
    
    v <- 1:n_row
    a<- 1 #amplitude in this case
    b<- 0
    T<- n_row
    phase <- 0
    use_cos <- FALSE
    vx <- sine_structure_fun(v,T,phase,a,b,use_cos)
    vx_rep <- unlist((lapply(1:n_row,FUN=function(j){rep(vx[j],time=n_col)})))  
    
    r2 <-setValues(r_init,vx_rep)  
    #plot(r2)
    
    r <- (r1+r2)*0.5  #combine patterns to make a elliptic surface, -0.5 could be modified
    #plot(r)
  }
  
  #change here to distance from edges..
  if(method=="use_edge"){ #does not work with large images
    #change...use gdal
    n_col <- ncol(r_init)
    n_row <- nrow(r_init)
    
    #xfrom
    r_init[1,1:n_col] <- 1
    r_init[n_row,1:n_col] <- 1
    r_init[1:n_row,1] <- 1
    r_init[1:n_row,n_col] <- 1
    #r_dist <- distance(r_init)
    #out_suffix_str
    
    ### Save raster file using tile id and Process ID
    srcfile <- file.path(out_dir_str,
                         paste("feature_target_",tile_no,"_",Sys.getpid(),out_suffix_str,file_format,sep=""))
    
    writeRaster(r_init,filename=srcfile,overwrite=T)
    #Sys.getpid
    dstfile <- file.path(out_dir_str,paste("feature_target_edge_distance",tile_no,"_",Sys.getpid(),out_suffix_str,file_format,sep=""))
    n_values <- "1"
    
    ### Generate distance from edge:
    cmd_str <- paste("gdal_proximity.py", 
                     srcfile, 
                     dstfile,
                     "-values",n_values,sep=" ")
    system(cmd_str)
    
    r_dist<- raster(dstfile) #read in the image created
    min_val <- cellStats(r_dist,min) 
    max_val <- cellStats(r_dist,max)
    r <- abs(r_dist - min_val)/ (max_val - min_val) #no need to inverse...
  } #too slow with R so used http://www.gdal.org/gdal_proximity.html
  
  if(method=="use_linear_weights"){
    #
    r_dist <- distance(r_init)
    min_val <- cellStats(r_dist,min) 
    max_val <- cellStats(r_dist,max)
    r <- abs(r_dist - max_val)/ (max_val - min_val)
  }
  
  ####### Now save layers for weights and prod weights
  #browser()
  extension_str <- extension(lf[i])
  raster_name_tmp <- gsub(extension_str,"",basename(lf[i]))
  raster_name <- file.path(out_dir_str,paste(raster_name_tmp,"_",method,"_weights_",out_suffix_str,file_format,sep=""))
  writeRaster(r, #write out the distance image
              NAflag=NA_flag_val,
              filename=raster_name,
              overwrite=TRUE)  
  
  ### Product weights computation
  r_var_prod <- r_in*r #this may be slow? could you GDAL instead
  raster_name_prod <- file.path(out_dir_str, paste(raster_name_tmp,"_",method,"_prod_weights_",out_suffix_str,file_format,sep=""))
  writeRaster(r_var_prod, 
              NAflag=NA_flag_val,
              filename=raster_name_prod,
              overwrite=TRUE)  
  
  weights_obj <- list(raster_name,raster_name_prod)
  names(weights_obj) <- c("r_weights","r_weights_prod")
  return(weights_obj)
}

mosaic_m_raster_list<-function(j,list_param){
  #This functions returns a subset of tiles from the modis grid.
  #Arguments: modies grid tile,list of tiles
  #Output: spatial grid data frame of the subset of tiles
  #Note that rasters are assumed to be in the same projection system!!
  #modified for global mosaic...still not working right now...
  
  #rast_list<-vector("list",length(mosaic_list))
  #for (i in 1:length(mosaic_list)){  
  # read the individual rasters into a list of RasterLayer objects
  # this may be changed so that it is not read in the memory!!!
  
  #parse output...
  
  #j<-list_param$j
  mosaic_list<-list_param$mosaic_list
  out_path<-list_param$out_path
  out_names<-list_param$out_rastnames
  file_format <- list_param$file_format
  NA_flag_val <- list_param$NA_flag_val
  out_suffix <- list_param$out_suffix
  ## Start
  
  if(class(mosaic_list[[j]])=="list"){
    m_list <- unlist(mosaic_list[[j]])
  }else{
    m_list <- mosaic_list[[j]]
  }
  input.rasters <- lapply(m_list, raster) #create raster image for each element of the list
  #inMemory(input.rasters[[1]])
  #note that input.rasters are not stored in memory!!
  mosaiced_rast<-input.rasters[[1]]
  
  for (k in 2:length(input.rasters)){
    mosaiced_rast<-mosaic(mosaiced_rast,input.rasters[[k]], tolerance=1,fun=mean)
    #mosaiced_rast<-mosaic(mosaiced_rast,raster(input.rasters[[k]]), fun=mean)
  }
  
  data_name<-paste("mosaiced_",sep="") #can add more later...
  #raster_name<-paste(data_name,out_names[j],".tif", sep="")
  raster_name<-paste(data_name,out_names[j],file_format, sep="")
  
  writeRaster(mosaiced_rast, NAflag=NA_flag_val,filename=file.path(out_path,raster_name),overwrite=TRUE)  
  #Writing the data in a raster file format...  
  rast_list<-file.path(out_path,raster_name)
  
  ## The Raster and rgdal packages write temporary files on the disk when memory is an issue. This can potential build up
  
  removeTmpFiles(h=0) #did not work if h is not set to 0
  ## end of remove section
  
  return(rast_list)
}

raster_match <- function(i,list_param){
  ### Read in parameters/arguments
  lf_files <- list_param$lf_files
  rast_ref <- list_param$rast_ref #name of reference file
  file_format <- list_param$file_format #".tif",".rst" or others
  python_bin <- list_param$python_bin
  out_suffix <- list_param$out_suffix
  out_dir_str <- list_param$out_dir_str
  
  ### START SCRIPT ##
  
  r_m <- raster(rast_ref) #ref image with resolution and extent to match
  
  set1f <- function(x){rep(NA, x)}
  
  inFilename <- lf_files[i]
  
  extension_str <- extension(inFilename)
  raster_name_tmp <- gsub(extension_str,"",basename(inFilename))
  #outFilename <- file.path(out_dir,paste(raster_name_tmp,"_","m_",out_suffix,file_format,sep="")) #for use in function later...
  
  raster_name <- file.path(out_dir_str,paste(raster_name_tmp,"_","m_",out_suffix,file_format,sep=""))#output file
  r_ref <- init(r_m, fun=set1f, filename=raster_name, overwrite=TRUE)
  #NAvalue(r_ref) <- -9999
  python_cmd <- file.path(python_bin,"gdalwarp")
  
  #cmd_str <- paste("/usr/bin/gdalwarp",inFilename,raster_name,sep=" ") #this may be a problem
  cmd_str <- paste(python_cmd,inFilename,raster_name,sep=" ") #this may be a problem
  #gdalwarp -t_srs '+proj=utm +zone=11 +datum=WGS84' raw_spot.tif utm11.tif
  system(cmd_str)
  
  ##return name of file created
  return(raster_name)
}

mosaicFiles <- function(lf_mosaic,mosaic_method="unweighted",num_cores=1,r_mask_raster_name=NULL,python_bin=NULL,mosaic_python="/nobackupp6/aguzman4/climateLayers/sharedCode/gdal_merge_sum_noDataTest.py",algorithm="R",match_extent=TRUE,df_points=NULL,NA_flag_val=-9999,file_format=".tif",out_suffix=NULL,out_dir=NULL,tmp_files=FALSE,data_type="Float32",scaling=NULL,values_range=NULL){
  #This functions mosaics tiles/files give a list of files. 
  #There are four options to mosaic:   use_sine_weights,use_edge,use_linear_weights, unweighted
  #Sine weights fits sine fuctions across rows and column producing elliptical/spherical patterns from center
  #Use edge uses the distance from the edge of the tiles/fies, higher weights towards the center
  #Linear weights use simple linear average from distance point feature (usually centroid)
  #Unweighted: average without and weigthing surface
  #In addition, option is given to the user to use R raster mosaic function or a python/gdal modified gdalmerge in mosaicing.
  #
  #INPUT Arguments: 
  #1)lf_mosaic: list of files to mosaic
  #2)mosaic_method: mosaic methods availbable:use_sine_weights,use_edge,use_linear_weights
  #3)num_cores: number of cores used in parallilization in mclapply
  #4)r_mask_raster_name: mask rference raster image
  #5)python_bin: location of python executables, defaut is NULL
  #6)mosaic_python: location/directory of python excecutable used for mosaicing with option sum/mean from Alberto Guzmann
  #7)df_points: point location used in weighting, defaul is NULL
  #8)NA_flag_val: raster flag values, e.g. -9999
  #9)file_format: raster format used, default is ".tif"
  #10)out_suffix: output suffix, default is NULL, it is recommended to add the variable name
  #             here e.g. dailyTmax and date!!
  #11)out_dir: output directory, default is NULL
  #12)algorithm: use R or python function
  #13)match extent: if TRUE match extent before mosaicing
  #14)tmp_files: if TRUE then keep temporary files
  #15)data_type: Float32 is default values for mosaicing
  #16)scaling: if NULL, use scaling 1, numeric value to multiply the values before conversion to integer
  #17)values_range: if NULL, don't screen
  #
  #OUTPUT:
  # Object is produced with 3 components:
  # 1) mean_mosaic: list of raster files from produced mosaic ,
  # 2) r_weights: list of raster files from weights 
  # 3) r_weights_prod: list of raster files from product weights (weights*value)
  # 4) method: weighting average used
  #
  
  ####################################
  ### BEGIN ####
  
  out_dir_str <- out_dir
  
  #if(tmp_files==T){
  out_suffix_str_tmp <- paste0(out_suffix,"_tmp")
  #}
  
  ### Maybe look at the data type format of input layers to mosaic
  #if(data_type==NULL){
  #  data_type <- "Int16" #should be a parameter!!
  #}else{
  #  data_type <- "Float32"
  #}
  if(is.null(scaling)){
    scaling <- 1
  }
  valid_range <- values_range #if NULL don't screen values!!
  #valid_range <- c(-100,100) #pass this as parameter!! (in the next update)
  if(data_type=="Int16"){
    data_type_str <- "INT2S"
  }
  lf_r_weights <- vector("list",length=length(lf_mosaic))
  
  rasterOptions(tmpdir=out_dir) #trying to control temporary files  written by the raster package

  #lf <- unlist(lf)
  
  ###############
  ### PART 2: prepare weights using tile rasters ############
  #methods availbable:use_sine_weights,use_edge,use_linear_weights
  
  if(mosaic_method=="use_linear_weights"){
    method <- "use_linear_weights"
    df_points <- df_centroids
    #df_points <- NULL
    r_feature <- NULL
    
    #lf <- list_param$lf
    #df_points <- list_param$df_points
    #r_feature <- list_param$r_feature #this should be change to a list
    #padding <- TRUE #if padding true then make buffer around edges??
    #method <- list_param$method #differnt methods available to create weights
    #NAflag,file_format,out_suffix etc...
    #NA_flag_val <- list_param$NA_flag
    #file_format <- list_param$file_format
    #out_suffix_str <- list_param$out_suffix_str
    #out_dir_str <- list_param$out_dir_str
    list_param_create_weights <- list(lf_mosaic,df_points,r_feature,method,NA_flag_val,file_format,out_suffix_str_tmp,out_dir_str) 
    names(list_param_create_weights) <- c("lf","df_points","r_feature","method","NA_flag","file_format","out_suffix_str","out_dir_str") 
    #num_cores <- 11
    
    #debug(create_weights_fun)
    #weights_obj <- create_weights_fun(1,list_param=list_param_create_weights)
    
    #This is the function creating the weights by tile. Distance from the centroids needs to be change from distance to
    #the edges...can use rows and columsn to set edges to 1 and 0 for the others.
    linear_weights_obj_list <- mclapply(1:length(lf_mosaic),FUN=create_weights_fun,list_param=list_param_create_weights,mc.preschedule=FALSE,mc.cores = num_cores)                           
    
    list_linear_r_weights <- lapply(1:length(linear_weights_obj_list), FUN=function(i,x){x[[i]]$r_weights},x=linear_weights_obj_list)
    list_linear_r_weights_prod <- lapply(1:length(linear_weights_obj_list), FUN=function(i,x){x[[i]]$r_weights_prod},x=linear_weights_obj_list)
    
    list_weights <- list_linear_r_weights
    list_weights_prod <- list_linear_r_weights_prod 
    
  }
  if(mosaic_method=="use_sine_weights"){
    
    ### Third use sine weights
    method <- "use_sine_weights"
    #df_points <- df_centroids
    df_points <- NULL
    r_feature <- NULL
    
    #list_param_create_weights <- list(lf_mosaic,df_points,r_feature,method,out_dir_str) 
    #names(list_param_create_weights) <- c("lf","df_points","r_feature","method","out_dir_str") 
    list_param_create_weights <- list(lf_mosaic,df_points,r_feature,method,NA_flag_val,file_format,out_suffix_str_tmp,out_dir_str) 
    names(list_param_create_weights) <- c("lf","df_points","r_feature","method","NA_flag","file_format","out_suffix_str","out_dir_str") 
    
    #num_cores <- 11
    
    #debug(create_weights_fun)
    #weights_obj <- create_weights_fun(1,list_param=list_param_create_weights)
    
    #This is the function creating the weights by tile. Distance from the centroids needs to be change from distance to
    #the edges...can use rows and columsn to set edges to 1 and 0 for the others.
    sine_weights_obj_list <- mclapply(1:length(lf_mosaic),FUN=create_weights_fun,list_param=list_param_create_weights,mc.preschedule=FALSE,mc.cores = num_cores)                           
    
    list_sine_r_weights <- lapply(1:length(sine_weights_obj_list), FUN=function(i,x){x[[i]]$r_weights},x=sine_weights_obj_list)
    list_sine_r_weights_prod <- lapply(1:length(sine_weights_obj_list), FUN=function(i,x){x[[i]]$r_weights_prod},x=sine_weights_obj_list)
    
    list_weights <- list_sine_r_weights
    list_weights_prod <- list_sine_r_weights_prod 
    
  }
  
  browser()
  if(mosaic_method=="use_edge_weights"){

    method <- "use_edge"
    df_points <- NULL
    r_feature <- NULL
    
    list_param_create_weights <- list(lf_mosaic,df_points,r_feature,method,NA_flag_val,file_format,out_suffix_str_tmp,out_dir_str) 
    names(list_param_create_weights) <- c("lf","df_points","r_feature","method","NA_flag","file_format","out_suffix_str","out_dir_str") 
    #list_param_create_weights <- list(lf_mosaic,df_points,r_feature,method,out_dir_str) 
    #names(list_param_create_weights) <- c("lf","df_points","r_feature","method","out_dir_str") 
    #num_cores <- 11
    #undebug(create_weights_fun)
    #weights_obj <- create_weights_fun(2,list_param=list_param_create_weights)
    
    #This is the function creating the weights by tile. Distance from the centroids needs to be change from distance to
    #the edges...can use rows and columsn to set edges to 1 and 0 for the others.
    num_cores <- 1
    
    use_edge_weights_obj_list <- mclapply(1:length(lf_mosaic),
                                          FUN=create_weights_fun,
                                          list_param=list_param_create_weights,
                                          mc.preschedule=FALSE,
                                          mc.cores = num_cores)
    #use_edge_weights_obj_list <- mclapply(2:length(lf_mosaic),
    #                                      FUN=create_weights_fun,
    #                                      list_param=list_param_create_weights,
    #                                      mc.preschedule=FALSE,mc.cores = num_cores) 
    
    #extract the list of files for weights and product weights
    list_edge_r_weights <- lapply(1:length(use_edge_weights_obj_list), FUN=function(i,x){x[[i]]$r_weights},x=use_edge_weights_obj_list)
    list_edge_r_weights_prod <- lapply(1:length(use_edge_weights_obj_list), FUN=function(i,x){x[[i]]$r_weights_prod},x=use_edge_weights_obj_list)
    
    #simplifly later...
    list_weights <- list_edge_r_weights
    list_weights_prod <- list_edge_r_weights_prod 
    #r_test <- raster(list_edge_r_weights[[1]])
    
  }
  
  ###############
  ### PART 3: prepare weightsfor mosaicing by matching extent ############
  
  ## Rasters tiles vary slightly in resolution, they need to be matched for the mosaic. Resolve issue in the 
  #mosaic function using gdal_merge to compute a reference image to mach.
  #This step of creating a merged raster can be avoided if a reference maks image is given
  #this needs to be changed to avoid further bugs!!!
  
  if(!is.null(r_mask_raster_name)){
    rast_ref <- r_mask_raster_name #the mask file is used as ref.
    #mask(raster(r_m_weighted_mean_raster_name),mask=r_mask_raster_name,filename=r_m_weighted_mean_mask_raster_name)
    #raster_name <- r_m_weighted_mean_mask_raster_name
  }else{
    rast_ref <- file.path(out_dir,paste("avg_",out_suffix,file_format,sep="")) #this is a the ref
    if(is.null(python_bin)){
      python_bin=""
    }
    
    python_cmd <- file.path(python_bin,"gdal_merge.py")  
    cmd_str <- paste("python",python_cmd,"-o ",rast_ref,paste(lf_mosaic,collapse=" ")) 
    system(cmd_str)
  }
  
  ## Create raster image for original predicted images with matching resolution and extent to the mosaic (reference image)
  
  #rast_ref <- file.path(out_dir,"avg.tif")
  r_ref <- raster(rast_ref)
  #plot(r_ref)
  
  if(mosaic_method%in%c("use_linear_weights","use_sine_weights","use_edge_weights")){
    
    if(algorithm=="python"){
      
      if(match_extent==TRUE){
        
        #If using R, it is necessary to match extent firt
        lf_files <- unlist(list_weights)
        
        ##Maching resolution is probably only necessary for the r mosaic function
        #Modify later to take into account option R or python...
        list_param_raster_match <- list(lf_files,rast_ref,file_format,python_bin,out_suffix_str_tmp,out_dir_str)
        names(list_param_raster_match) <- c("lf_files","rast_ref","file_format","python_bin","out_suffix","out_dir_str")
        
        #undebug(raster_match)
        #r_test <- raster_match(1,list_param_raster_match)
        #r_test <- raster(raster_match(1,list_param_raster_match))
        
        list_weights_m <- mclapply(1:length(lf_files),FUN=raster_match,list_param=list_param_raster_match,mc.preschedule=FALSE,mc.cores = num_cores)                           
        
        lf_files <- unlist(list_weights_prod)
        list_param_raster_match <- list(lf_files,rast_ref,file_format,python_bin,out_suffix_str_tmp,out_dir_str)
        names(list_param_raster_match) <- c("lf_files","rast_ref","file_format","python_bin","out_suffix","out_dir_str")
        
        #num_cores <-11
        list_weights_prod_m <- mclapply(1:length(lf_files),FUN=raster_match,list_param=list_param_raster_match,mc.preschedule=FALSE,mc.cores = num_cores)                           
        
      }else{
        list_weights_m <- list_weights
        list_weights_prod_m <- list_weights_prod 
      }
      
      #pattern_str <- paste("*.","predicted_mod1",".*.",day_to_mosaic[i],".*.tif",sep="")
      #lf_day_to_mosaic <- lapply(1:length(unlist(in_dir_mosaics)),FUN=function(k){list.files(path=unlist(in_dir_mosaics)[k],pattern=pattern_str,full.names=T,recursive=T)}) 
      #lf_day_to_mosaic <- unlist(lf_day_to_mosaic)
      #write.table(lf_day_to_mosaic,file=file.path(out_dir,paste("list_to_mosaics_",day_to_mosaic[i],".txt",sep="")))
      #filename_list_mosaics <- file.path(out_dir,paste("list_to_mosaics_",day_to_mosaic[i],".txt",sep=""))
      
      filename_list_mosaics_weights_m <- file.path(out_dir_str,
                                                   paste("list_to_mosaics_","weights_",mosaic_method,"_",out_suffix_str_tmp,".txt",sep=""))
      filename_list_mosaics_prod_weights_m <- file.path(out_dir_str,
                                                        paste("list_to_mosaics_","prod_weights_",mosaic_method,"_",out_suffix_str_tmp,".txt",sep=""))
      
      #writeLines(unlist(list_weights_m),con=filename_list_mosaics_weights_m) #weights files to mosaic 
      #writeLines(unlist(list_weights_prod_m),con=filename_list_mosaics_prod_weights_m) #prod weights files to mosaic
      
      #writeLines(unlist(list_weights_m),
      #           con=filename_list_mosaics_weights_m) #weights files to mosaic 
      write.table(as.data.frame(unlist(list_weights_m)),
                   file=filename_list_mosaics_weights_m,
                   sep=" ",
                   col.names=F,
                  row.names=F)
      
      #writeLines(unlist(list_weights_prod_m),con=filename_list_mosaics_prod_weights_m) #prod weights files to mosaic
      
      write.table(as.data.frame(unlist(list_weights_prod_m)),
                  file=filename_list_mosaics_prod_weights_m,
                  sep=" ",
                  col.names=F,
                  row.names=F)
      
      #out_mosaic_name_weights_m <- r_weights_sum_raster_name <- file.path(out_dir,paste("r_weights_sum_m_",mosaic_method,"_weighted_mean_",out_suffix,".tif",sep=""))
      #out_mosaic_name_prod_weights_m <- r_weights_sum_raster_name <- file.path(out_dir,paste("r_prod_weights_sum_m_",mosaic_method,"_weighted_mean_",out_suffix,".tif",sep=""))
      out_mosaic_name_weights_m  <- file.path(out_dir_str,paste("r_weights_sum_m_",mosaic_method,"_weighted_mean_",out_suffix_str_tmp,".tif",sep=""))
      out_mosaic_name_prod_weights_m <- file.path(out_dir_str,paste("r_prod_weights_sum_m_",mosaic_method,"_weighted_mean_",out_suffix_str_tmp,".tif",sep=""))
      
      module_path <- mosaic_python #this should be a parameter for the function...
      #browser()
      #python /nobackupp6/aguzman4/climateLayers/sharedCode/gdal_merge_sum.py 
      #--config GDAL_CACHEMAX=1500 --overwrite=TRUE -o  outputname.tif --optfile input.txt
      if(!is.null(r_mask_raster_name)){
        #d
        #different extent between mask and output if match extent is false!!
        #match resolution and extent first
        
        #lf_files <- c(r_m_weighted_mean_raster_name) #file(s) to be matched
        #lf_files <- c(r_m_weighted_mean_raster_name) #match to mask
        rast_ref_name <- r_mask_raster_name
        #r_mask <- raster(rast_ref)
        #extent_r_mask <- as.numeric(as.matrix(extent(r_mask)))
        #res_pix <- res(r_mask)
        #c(xmin,ymax,xmax,ymin)
        #c(ulx,uly,lrx,lry)
        #extent_str <- c(extent_r_mask[1],extent_r_mask[4],extent_r_mask[3],extent_r_mask[2])
      }else{
        rast_ref_name <- NULL
      }
      #cmd_mosaic_gam_CAI_dailyTmax_19840101_reg1_1984.txt
      #python /nobackupp6/aguzman4/climateLayers/sharedCode//gdal_merge_sum.py --config GDAL_CACHEMAX=1500 --overwrite=TRUE -o /nobackupp8/bparmen1/climateLayers/out
      #debug(mosaic_python_merge)
      mosaic_weights_obj <- mosaic_python_merge(NA_flag_val=NA_flag_val,
                                                module_path=mosaic_python,
                                                module_name="gdal_merge_sum.py",
                                                input_file=filename_list_mosaics_weights_m,
                                                out_mosaic_name=out_mosaic_name_weights_m,
                                                raster_ref_name = rast_ref_name) ##if NA, not take into account
      r_weights_sum_raster_name <- mosaic_weights_obj$out_mosaic_name
      cmd_str1 <- mosaic_weights_obj$cmd_str
      #r_prod_sum_raster_name <- mosaic_python_merge(module_path=mosaic_python,
      #                                              module_name="gdal_merge_sum.py",
      #                                              input_file=filename_list_mosaics_prod_weights_m,
      #                                              out_mosaic_name=out_mosaic_name_prod_weights_m)
      
      mosaic_prod_weights_obj <- mosaic_python_merge(NA_flag_val=NA_flag_val,
                                                     module_path=mosaic_python,
                                                     module_name="gdal_merge_sum.py",
                                                     input_file=filename_list_mosaics_prod_weights_m,
                                                     out_mosaic_name=out_mosaic_name_prod_weights_m,
                                                     raster_ref_name = rast_ref_name)
      r_prod_sum_raster_name <- mosaic_prod_weights_obj$out_mosaic_name
      cmd_str2 <- mosaic_prod_weights_obj$cmd_str
      #write out python command used for mosaicing
      cmd_mosaic_logfile <- file.path(out_dir,paste("cmd_mosaic_",out_suffix,".txt",sep=""))
      writeLines(cmd_str1,con=cmd_mosaic_logfile) #weights files to mosaic 
      #writeLines(cmd_str2,con=file.path(out_dir,paste("cmd_mosaic_",out_suffix,".txt",sep=""))) #weights files to mosaic 
      cat(cmd_str2, file=cmd_mosaic_logfile, append=TRUE, sep = "\n")
    }
    
    if(algorithm=="R"){
      
      #If using R, it is necessary to match extent firt
      
      lf_files <- unlist(list_weights)
      
      ##Maching resolution is probably only necessary for the r mosaic function
      #MOdify later to take into account option R or python...
      list_param_raster_match <- list(lf_files,rast_ref,file_format,python_bin,out_suffix_str_tmp,out_dir_str)
      names(list_param_raster_match) <- c("lf_files","rast_ref","file_format","python_bin","out_suffix","out_dir_str")
      
      #undebug(raster_match)
      #r_test <- raster(raster_match(1,list_param_raster_match))
      
      list_weights_m <- mclapply(1:length(lf_files),FUN=raster_match,list_param=list_param_raster_match,mc.preschedule=FALSE,mc.cores = num_cores)                           
      
      lf_files <- unlist(list_weights_prod)
      list_param_raster_match <- list(lf_files,rast_ref,file_format,python_bin,out_suffix_str_tmp,out_dir_str)
      names(list_param_raster_match) <- c("lf_files","rast_ref","file_format","python_bin","out_suffix","out_dir_str")
      
      #num_cores <-11
      list_weights_prod_m <- mclapply(1:length(lf_files),FUN=raster_match,list_param=list_param_raster_match,mc.preschedule=FALSE,mc.cores = num_cores)                           
      
      
      #The file to do the merge is /nobackupp6/aguzman4/climateLayers/sharedCode/gdal_merge_sum.py. Sample call below.
      #python /nobackupp6/aguzman4/climateLayers/sharedCode/gdal_merge_sum.py --config GDAL_CACHEMAX=1500 --overwrite=TRUE -o  outputname.tif --optfile input.txt
      ##Make this a function later
      #list_weights_m <- list(list_linear_weights_m,list_edge_weights_m,list_sine_weights_m)
      #list_weights_prod_m <- list(list_linear_weights_prod_m,list_edge_weights_prod_m,list_sine_weights_prod_m)
      #list_methods <- c("linear","edge","sine")
      list_mosaiced_files <- vector("list",length=1)
      
      list_args_weights <- list_weights_m
      list_args_weights_prod <- list_weights_prod_m
      method_str <- method
      
      #making a list of raster object before mosaicing
      list_args_weights <- lapply(1:length(list_args_weights), FUN=function(i,x){raster(x[[i]])},x=list_args_weights)
      
      #get the list of weights product into raster objects
      
      list_args_weights_prod <- lapply(1:length(list_args_weights_prod), FUN=function(i,x){raster(x[[i]])},x=list_args_weights_prod)
      list_args_weights_prod$fun <- "sum" #use sum while mosaicing
      list_args_weights_prod$na.rm <- TRUE #deal with NA by removal
      r_weights_sum_raster_name <- file.path(out_dir_str,paste("r_weights_sum_m_",method_str,"_weighted_mean_",out_suffix_str_tmp,".tif",sep=""))
      list_args_weights$filename <- r_weights_sum_raster_name
      list_args_weights$overwrite<- TRUE
      list_args_weights_prod$overwrite<- TRUE  #add to overwrite existing image  
      
      list_args_weights$fun <- "sum" #we want the sum to compute the weighted mean
      list_args_weights$na.rm <- TRUE
      r_prod_sum_raster_name <- file.path(out_dir_str,paste("r_prod_sum_m_",method_str,"_weighted_mean_",out_suffix_str_tmp,".tif",sep=""))
      list_args_weights_prod$filename <- r_prod_sum_raster_name
      
      #Mosaic files: this is where we can use Alberto Python function but modified with option for
      #sum in addition ot the current option for mean.
      #This took 23 minutes!
      r_weights_sum <- do.call(mosaic,list_args_weights) #weights sum image mosaiced
      #This took 23 minutes with the R function
      r_prod_sum <- do.call(mosaic,list_args_weights_prod) #weights sum product image mosacied
      
    }
    
    
    #r_m_weighted_mean <- r_prod_sum/r_weights_sum #this is the mosaic using weighted mean...
    
    r_m_weighted_mean_raster_name <- file.path(out_dir_str,paste("r_m_",mosaic_method,"_weighted_mean_",out_suffix,".tif",sep=""))
    #r_m_weighted_mean_raster_name <- "test_tmp.tif"
    if(is.null(python_bin)){
      python_bin=""
    }
    
    ##Add here the int and scaling?
    #note that the nodata was fixed...
    #if not null use the value specificied in the parameters
    #browser()
    python_cmd <- file.path(python_bin,"gdal_calc.py")
    cmd_str3 <- paste(python_cmd, 
                      paste("-A ", r_prod_sum_raster_name,sep=""),
                      paste("-B ", r_weights_sum_raster_name,sep=""),
                      paste("--outfile=",r_m_weighted_mean_raster_name,sep=""),
                      paste("--type=",data_type,sep=""),
                      "--co='COMPRESS=LZW'",
                      paste("--NoDataValue=",NA_flag_val,sep=""),
                      paste("--calc='(A/B)*",scaling,"'",sep=""),
                      "--overwrite",sep=" ") #division by zero is problematic...
    system(cmd_str3)
    
    ## Skipping this step now...
    #if(!is.null(r_mask_raster_name)){
    #different extent between mask and output if match extent is false!!
    #match resolution and extent first
    
    #   #lf_files <- c(r_m_weighted_mean_raster_name) #file(s) to be matched
    #  lf_files <- c(r_m_weighted_mean_raster_name) #match to mask
    #    rast_ref <- r_mask_raster_name
    #    raster(rast_ref)
    #    extent_r_mask <- extent(r_mask)
    ##Maching resolution is probably only necessary for the r mosaic function
    #Modify later to take into account option R or python...
    #   list_param_raster_match <- list(lf_files,rast_ref,file_format,python_bin,out_suffix,out_dir)
    #    names(list_param_raster_match) <- c("lf_files","rast_ref","file_format","python_bin","out_suffix","out_dir_str")
    
    #undebug(raster_match)
    #   r_m_weighted_mean_raster_name_matched <- raster_match(1,list_param_raster_match)
    #}
    
    #writeRaster(r_m_weighted_mean, NAflag=NA_flag_val,filename=raster_name,overwrite=TRUE)  
    
    ###Starting rescaling, switched by masking first then screening to avoid potential problems
    ##Can merge one and 2 with parentheses operations!!, make this a function?
    ##Reclassify with valid range: -100,100
    
    raster_name <- r_m_weighted_mean_raster_name
    #raster_name <- r_m_weighted_mean_raster_name_matched
    max_val <- valid_range[2]*scaling #set min_valid
    raster_name_rec1 <- file.path(out_dir_str,paste("r_m_",mosaic_method,"_weighted_mean_rec1_",out_suffix,"_tmp",".tif",sep=""))
    #rec_tmp1 <- file.path(out_dir_str,paste("r_m_",mosaic_method,"_weighted_mean_rec_",out_suffix,".tif",sep=""))
    cmd_str4 <- paste(python_cmd, 
                      paste("-A ", raster_name,sep=""),
                      paste("--outfile=",raster_name_rec1,sep=""),
                      paste("--type=",data_type,sep=""),
                      paste("--NoDataValue=",NA_flag_val,sep=""),
                      "--co='COMPRESS=LZW'",
                      paste("--calc='(A>",max_val,")*",NA_flag_val,"'",sep=""),
                      "--overwrite",sep=" ") #division by zero is problematic...
    system(cmd_str4)
    
    raster_name_rec2 <- file.path(out_dir_str,paste("r_m_",mosaic_method,"_weighted_mean_rec2_",out_suffix,"_tmp",".tif",sep=""))
    min_val <- valid_range[1]*scaling #set min_valid as a input
    cmd_str5 <- paste(python_cmd, 
                      paste("-A ", raster_name,sep=""),
                      paste("--outfile=",raster_name_rec2,sep=""),
                      paste("--type=",data_type,sep=""),
                      "--co='COMPRESS=LZW'",
                      paste("--NoDataValue=",NA_flag_val,sep=""),
                      paste("--calc='(A<",min_val,")*",NA_flag_val,"'",sep=""),
                      "--overwrite",sep=" ") #division by zero is problematic...
    system(cmd_str5)
    
    ##Now combine A and B and multiply by mask??
    r_m_weighted_mean_raster_name_rec <- file.path(out_dir_str,paste("r_m_",mosaic_method,"_weighted_mean_masked_",out_suffix,".tif",sep=""))
    cmd_str6 <- paste(python_cmd, 
                      paste("-A ", raster_name_rec1,sep=""),
                      paste("-B ", raster_name_rec2,sep=""),
                      paste("-C ", raster_name,sep=""), #if mask exists
                      paste("--outfile=",r_m_weighted_mean_raster_name_rec,sep=""),
                      paste("--type=",data_type,sep=""),
                      "--co='COMPRESS=LZW'",
                      paste("--NoDataValue=",NA_flag_val,sep=""),
                      paste("--calc='(((((A+B))<",min_val,")*",NA_flag_val,")+1)*C'",sep=""),
                      "--overwrite",sep=" ") #division by zero is problematic...
    system(cmd_str6)    
    
    #check if file exists first...
    
    ### End of rescaling section
    #r_m_use_edge_weights_weighted_mean_rec_gam_CAI_dailyTmax_19910101_reg4_tmp.tif
    
    #browser() #22minutes for one mosaic
    
    cmd_mosaic_logfile <- file.path(out_dir,paste("cmd_mosaic_",out_suffix,".txt",sep=""))
    
    writeLines(cmd_str1,con=cmd_mosaic_logfile) #weights files to mosaic 
    #writeLines(cmd_str2,con=file.path(out_dir,paste("cmd_mosaic_",out_suffix,".txt",sep=""))) #weights files to mosaic 
    cat(cmd_str2, file=cmd_mosaic_logfile, append=TRUE, sep = "\n")
    cat(cmd_str3, file=cmd_mosaic_logfile, append=TRUE, sep = "\n")
    cat(cmd_str4, file=cmd_mosaic_logfile, append=TRUE, sep = "\n")
    cat(cmd_str5, file=cmd_mosaic_logfile, append=TRUE, sep = "\n")
    cat(cmd_str6, file=cmd_mosaic_logfile, append=TRUE, sep = "\n")
    
    #cmd_str <- "/nobackupp6/aguzman4/climateLayers/sharedModules/bin/gdal_calc.py -A r_prod_weights_sum_m_use_edge_weights_weighted_mean_gam_CAI_dailyTmax_19920101_reg4_run10_1500x4500_global_analyses_pred_1992_10052015.tif -B r_weights_sum_m_use_edge_weights_weighted_mean_gam_CAI_dailyTmax_19920101_reg4_run10_1500x4500_global_analyses_pred_1992_10052015.tif --outfile='test2.tif' --calc='A/B' --overwrite"
    
    #now use the mask
    if(!is.null(r_mask_raster_name)){
      #different extent between mask and output if match extent is false!!
      #match resolution and extent first
      
      #lf_files <- c(r_m_weighted_mean_raster_name) #file(s) to be matched
      #lf_files <- c(r_m_weighted_mean_raster_name_rec) #match to mask
      rast_ref <- r_mask_raster_name
      ##Maching resolution is probably only necessary for the r mosaic function
      #Modify later to take into account option R or python...
      #list_param_raster_match <- list(lf_files,rast_ref,file_format,python_bin,out_suffix,out_dir)
      #names(list_param_raster_match) <- c("lf_files","rast_ref","file_format","python_bin","out_suffix","out_dir_str")
      
      #undebug(raster_match)
      #r_m_weighted_mean_raster_name_matched <- raster_match(1,list_param_raster_match)
      #output
      r_m_weighted_mean_mask_raster_name <- file.path(out_dir_str,paste("r_m_",mosaic_method,"_weighted_mean_mask_",out_suffix,".tif",sep=""))
      
      in_raster_name <- r_m_weighted_mean_raster_name_rec
      mask(raster(in_raster_name),
           mask=raster(r_mask_raster_name),
           filename=r_m_weighted_mean_mask_raster_name,
           datatype=data_type_str,
           options=c("COMPRESS=LZW"),
           overwrite=TRUE,
           NAflag=NA_flag_val)
      
      ##Now combine A and B and multiply by mask?? This must be faster than the mask option in R
      #The mask must be in 1,NA format with 1 being valid values being considered in roi.
      #r_m_weighted_mean_raster_name_rec <- file.path(out_dir_str,paste("r_m_",mosaic_method,"_weighted_mean_rec_",out_suffix,"_tmp",".tif",sep=""))
      #browser()
      #cmd_mosaic_logfile <- file.path(out_dir,paste("cmd_mosaic_",out_suffix,".txt",sep=""))
      #in_raster_name <- r_m_weighted_mean_raster_name_matched
      
      
      #cmd_str7 <- paste(python_cmd, 
      #               paste("-A ", in_raster_name,sep=""),
      #               paste("-B ", r_mask_raster_name,sep=""),
      #               paste("--outfile=",r_m_weighted_mean_mask_raster_name,sep=""),
      #               paste("--type=",data_type,sep=""),
      #               "--co='COMPRESS=LZW'",
      #               paste("--NoDataValue=",NA_flag_val,sep=""),
      #               paste("--calc='(A*B)'",sep=""),
      #               "--overwrite",sep=" ") #division by zero is problematic...
      #system(cmd_str7)  
      #cat(cmd_str7, file=cmd_mosaic_logfile, append=TRUE, sep = "\n")
      
      ##Set min max and NA value
      r_mosaiced <- raster(r_m_weighted_mean_mask_raster_name)
      r_mosaiced <- setMinMax(r_mosaiced)
      NAvalue(r_mosaiced) <- NA_flag_val
      raster_name <- r_m_weighted_mean_mask_raster_name
      #browser()
      #-32,768 is NA
    }else{
      raster_name <- r_m_weighted_mean_raster_name
    }
    
  } #end of weighted
  
  ###
  if(mosaic_method=="unweighted"){
    #### Fourth use original images
    #macth file to mosaic extent using the original predictions
    
    if(match_extent==TRUE){
      lf_files <- lf_mosaic
      list_param_raster_match <- list(lf_files,rast_ref,file_format,out_suffix,out_dir)
      names(list_param_raster_match) <- c("lf_files","rast_ref","file_format","out_suffix","out_dir_str")
      list_pred_m <- mclapply(1:length(lf_files),FUN=raster_match,list_param=list_param_raster_match,mc.preschedule=FALSE,mc.cores = num_cores)                           
    }else{
      list_pred_m <- lf_mosaic
    }
    #list_mosaiced_files <- list.files(pattern="r_m.*._weighted_mean_.*.tif")
    
    #names(list_mosaiced_files) <- c("edge","linear","sine")
    
    #### NOW unweighted mean mosaic
    
    #get the original predicted image to raster (matched previously to the mosaic extent)
    list_args_pred_m <- list_pred_m
    #list_args_pred_m <- (mixedsort(list.files(pattern="^gam_CAI.*.m_mosaic_run10_1500x4500_global_analyses_03252015.tif")))
    list_args_pred_m <- lapply(1:length(list_args_pred_m), FUN=function(i,x){raster(x[[i]])},x=list_args_pred_m)
    
    list_args_pred_m$fun <- "mean"
    list_args_pred_m$na.rm <- TRUE
    list_args_pred_m$overwrite<- TRUE  #add to overwrite existing image 
    #list_args_pred_m$filename <- 
    
    #Mosaic files using R raster mosaic 
    r_m_mean <- do.call(mosaic,list_args_pred_m) #this is unweighted mean from the predicted raster
    
    r_m_mean_raster_name <- file.path(out_dir,paste("r_m_mean_",out_suffix,".tif",sep=""))
    writeRaster(r_m_mean, NAflag=NA_flag_val,filename=r_m_mean_raster_name,overwrite=TRUE)  #unweighted mean
    
    #r_m_mean_unweighted <- paste("r_m_mean_",out_suffix,".tif",sep="")
    #list_weights <- NULL
    #list_weights_prod <- NULL
    
    
    if(!is.null(r_mask_raster_name)){
      #different extent between mask and output if match extent is false!!
      #match resolution and extent first
      
      lf_files <- c(r_m_mean_raster_name) #file(s) to be matched
      rast_ref <- r_mask_raster_name
      ##Maching resolution is probably only necessary for the r mosaic function
      #Modify later to take into account option R or python...
      list_param_raster_match <- list(lf_files,rast_ref,file_format,python_bin,out_suffix,out_dir)
      names(list_param_raster_match) <- c("lf_files","rast_ref","file_format","python_bin","out_suffix","out_dir_str")
      
      #undebug(raster_match)
      r_m_mean_raster_name_matched <- raster_match(1,list_param_raster_match)
      
      r_m_mean_mask_raster_name <- file.path(out_dir,paste("r_m_",method_str,"_unweighted_mean_mask_",out_suffix,".tif",sep=""))
      mask(raster( r_m_mean_raster_name_matched),mask=raster(r_mask_raster_name),
           filename=r_m_mean_mask_raster_name,overwrite=TRUE)
      raster_name <- r_m_mean_mask_raster_name
      r_raster <- raster(raster_name)
      r_raster <- setMinMax(r_raster) #set correct min and max in the file
    }else{
      raster_name <- r_m_mean_raster_name
      r_raster <- raster(raster_name)
      r_raster <- setMinMax(r_raster) #set correct min and max in the file
    }
    
  }
  
  ########## clean up the disk/directories before ending the function ####
  
  if(tmp_files==F){ #if false...delete all files with "_tmp"
    lf_tmp <- list.files(pattern="*.*tmp*.*",path=out_dir_str,full.names=T)
    ##now delete temporary files...
    file.remove(lf_tmp)
  }
  
  ### remove temporary files generated by the raster package:
  #raster_tmp_2017-02-19_134558_60103_03411.gri #this can be up to 1.4gb for reg1 or reg4 
  lf_raster_tmp <- list.files(pattern="^raster_tmp_*.*",path=out_dir_str,full.names=T)
  file.remove(lf_raster_tmp)
  
  #Create return object
  mosaic_obj <- list(raster_name,list_weights,list_weights_prod,mosaic_method)
  names(mosaic_obj) <- c("mean_mosaic","r_weights","r_weigths_prod","method")
  save(mosaic_obj,file=file.path(out_dir_str,paste(mosaic_method,"_","mosaic_obj_",out_suffix,".RData",sep="")))
  return(mosaic_obj)
}



remove_errors_list<-function(list_items){
  
  #This function removes "error" items in a list
  list_tmp<-list_items
  if(is.null(names(list_tmp))){
    names(list_tmp) <- paste("l",1:length(list_tmp),sep="_")
    names(list_items) <- paste("l",1:length(list_tmp),sep="_")
  }
  
  for(i in 1:length(list_items)){
    if(inherits(list_items[[i]],"try-error")){
      list_tmp[[i]]<-0
    }else{
      list_tmp[[i]]<-1
    }
  }
  cnames<-names(list_tmp[list_tmp>0])
  x <- list_items[match(cnames,names(list_items))]
  return(x)
}

get_mosaic_files_fun  <- function(i,day_to_mosaic,in_dir_tiles_tmp,year_processed){
  #Find relevant files to mosaic for world/global mosaicing (by regions)
  lf_tmp <- lapply(1:length(in_dir_tiles_tmp),
                   FUN=function(j){
                     searchStr = paste(in_dir_tiles_tmp[j],"/","output_*",year_processed,"/r_m_use_edge_weights_weighted_mean_mask_","*",day_to_mosaic[i],"*.tif",sep="")
                     Sys.glob(searchStr)})
  lf_tmp <- unlist(lf_tmp)
  return(lf_tmp)
}

##################### END OF SCRIPT ######################
