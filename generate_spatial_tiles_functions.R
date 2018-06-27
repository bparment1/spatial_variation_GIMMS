############### Spatial utility: General code for generating tiles/spatial subsets for AREA  ########## 
## 
## DATE CREATED: 06/08/2017
## DATE MODIFIED: 06/27/2018
## AUTHORS: Benoit Parmentier 
## Version: 1
## PROJECT: General purpose
## ISSUE: 
## TO DO: Make this a function later
##
## COMMIT: Adding overlapping in tiles generation
##
## Links to investigate:
#
###################################################
#

###### Library used

library(gtools)                              # loading some useful tools 
library(sp)                                  # Spatial pacakge with class definition by Bivand et al.
library(spdep)                               # Spatial pacakge with methods and spatial stat. by Bivand et al.
library(rgdal)                               # GDAL wrapper for R, spatial utilities
library(gdata)                               # various tools with xls reading, cbindX
library(rasterVis)                           # Raster plotting functions
library(parallel)                            # Parallelization of processes with multiple cores
library(maptools)                            # Tools and functions for sp and other spatial objects e.g. spCbind
library(maps)                                # Tools and data for spatial/geographic objects
library(plyr)                                # Various tools including rbind.fill
library(spgwr)                               # GWR method
library(rgeos)                               # Geometric, topologic library of functions
library(gridExtra)                           # Combining lattice plots
library(colorRamps)                          # Palette/color ramps for symbology
library(ggplot2)

###### Functions used in this script

#This function creates a spatial polygon data frame object for the extent matching a raster input
create_polygon_from_extent<-function(reg_ref_rast,outDir=NULL,outSuffix=NULL){
  #This functions returns polygon sp from input rast
  #Input Arguments: 
  #reg_ref_rast: input ref rast
  #outDir : output directory, if NULL then the current dir in used
  #outSuffix: output suffix used for the naming of the shapefile
  #Output: 
  #reg_outline_poly: spatial polygon data.frame
  #
  if(is.null(outDir)){
    outDir=getwd()
  }
  if(is.null(outSuffix)){
    outSuffix=""
  }
  ref_e <- extent(reg_ref_rast) #extract extent from raster object
  reg_outline_poly <- as(ref_e, "SpatialPolygons") #coerce raster extent object to SpatialPolygons from sp package 
  reg_outline_poly <- as(reg_outline_poly, "SpatialPolygonsDataFrame") #promote to spdf
  proj4string(reg_outline_poly) <- projection(reg_ref_rast) #Assign projection to spdf
  infile_reg_outline <- paste("reg_out_line_",out_suffix,".shp",sep="") #name of newly crated shapefile with the extent
  writeOGR(reg_outline_poly,dsn= outDir,layer= sub(".shp","",infile_reg_outline), 
           driver="ESRI Shapefile",overwrite_layer="TRUE")
  
  return(reg_outline_poly) #return spdf
}

generate_grid_tiles <- function(ref_file,n_tile,n_tile_x, n_tile_y,out_suffix=NULL,out_dir=NULL){
  
  r <- raster(ref_file)
  
  # for the time being generate a non-overlapping grid tiling and crop
  extent_val <- extent(r)
  bbox_val <- st_bbox(r)
  test_sp <- as(extent_val, 'SpatialPolygons')
  outline_sf <-as(test_sp,"sf")
  
  #Can buffer?
  
  #test_grid <- st_make_grid(outline_sf, n=18)
  tile_grid <- st_make_grid(outline_sf, n=n_tile)
  
  plot(r)
  plot(tile_grid,add=T)
  #plot(test_grid[56],add=T,col="red")
  
  if(!is.null(out_suffix)){
    out_suffix_s <- paste0("_",out_suffix)
  }else{
    out_suffix_s <- ""
  }
  out_filename <- paste0("grid_tiles",out_suffix_s,".shp")
  out_grid_filename <- file.path(out_dir,out_filename)
  st_write(tile_grid,dsn=out_grid_filename)
  
  #Generate overlapping grid option to come later
  
  return(tile_grid)
}


generate_tiles_from_extent <- function(r_ref,y_ratio=2,x_ratio=2,x_overlap=0,y_overlap=0,out_suffix=NULL,out_dir=NULL){
  #
  r_ref <- raster(infile_ref_r)
  plot(r_ref)
  
  reg_outline_sp <- create_polygon_from_extent(reg_ref_rast=r_ref,outDir=NULL,outSuffix=NULL)
  plot(reg_outline_sp,add=T,border="red")
  
  reg_centroid <- gCentroid(reg_outline_sp)
  
  reg_extent <- extent( reg_outline_sp) #get boudning box of extent
  
  #Now make this an extent and a polygon, then shift the polygon on the right by half
  
  #First create one tile:
  
  xmin_new <- xmin(reg_extent) #Upper left corner
  xmax_new <- xmin(reg_extent) + ((xmax(reg_extent)- xmin(reg_extent))/x_ratio) #upper right corner
  ymin_new <- ymax(reg_extent)-((ymax(reg_extent)- ymin(reg_extent))/y_ratio)
  ymax_new <- ymax(reg_extent)
  
  range_e_x <- ((xmax(reg_extent)- xmin(reg_extent))/x_ratio)
  range_e_y <- ((ymax(reg_extent)- ymin(reg_extent))/y_ratio)
  
  #xmin,xmax,ymin,ymax
  tile_1_e <- extent(xmin_new,xmax_new,ymin_new,ymax_new)
  tile_1_extent_sp <- as(tile_1_e, "SpatialPolygons") #coerce raster extent object to SpatialPolygons from sp package 
  tile_1_extent_spdf <- as(tile_1_extent_sp, "SpatialPolygonsDataFrame") #promote to spdf
  proj4string(tile_1_extent_sp) <- projection(r_ref) #Assign projection to spdf
  
  ###### Now generate additional tiles:
  
  
  ##### No overlap
  
  if(x_overlap==0 || y_overlap==0){
    n_tiles <- x_ratio*y_ratio 
    
    list_row <- vector("list",length=y_ratio)
    
    list_col <- vector("list",length=x_ratio)
    
    #j <- 1
    for(j in 1:y_ratio){
      
      for(i in 1:x_ratio){
        
        x_shift <- range_e_x*(i-1)
        y_shift <- -1*range_e_y*(j-1)
        
        tile_extent_spdf <- shift(tile_1_extent_spdf,
                                  x= x_shift,
                                  y=y_shift)
        list_col[[i]] <- tile_extent_spdf
      }
      list_row[[j]] <- list_col
    }
    
    list_matrix_tiles <- list_row
    list_tiles<- unlist(list_matrix_tiles)
    
  }
    
  #With overlap x and overlap y options later on
  
  if(x_overlap!=0 || y_overlap!=0){
    
    n_x <- round(x_ratio * ((1 - x_overlap)^-1))
    n_y <- round(y_ratio * ((1 - y_overlap)^-1))
    
    n_tiles <- n_x*n_y 
    
    list_row <- vector("list",length=n_y)
    
    list_col <- vector("list",length=n_x)
    
    #j <- 1
    for(j in 1:n_y){
      
      for(i in 1:n_x){
        
        x_shift <- range_e_x*(i-1) * (1-x_overlap)
        y_shift <- -1*range_e_y*(j-1) * (1-y_overlap)
        
        tile_extent_spdf <- shift(tile_1_extent_spdf,
                                  x= x_shift,
                                  y=y_shift)
        list_col[[i]] <- tile_extent_spdf
      }
      list_row[[j]] <- list_col
    }
    
    list_matrix_tiles <- list_row
    list_tiles<- unlist(list_matrix_tiles)
    
  }
    
  ### Can clean up list of tiles based on overlap with extent afterwards!! It could be based on minimum overlap with the 
  ## extentent area
  
  plot(r_ref)
  plot(reg_outline_sp,add=T,border="red")
  
  list_centroids <- vector("list",length=n_tiles)
  for(k in 1:n_tiles){
    tile_spdf <- list_tiles[[k]]
    centroid_sp <- gCentroid(tile_spdf)
    
    plot(tile_spdf,border="black",add=T)
    text(coordinates(centroid_sp)[1],coordinates(centroid_sp)[2],k)
    list_centroids[[k]] <- centroid_sp
  }
  
  out_suffix_s<- out_suffix
  list_tiles_names <- paste("tile_",1:length(list_tiles),"_",out_suffix,".shp",sep="")
  
  ### create a dir to store the shapefiles
  if(!is.null(out_dir)){
    out_dir_s <- out_dir
  }else{
    out_dir_s <- "tiles"
    if(!dir.exists(out_dir_s)){
      dir.create(out_dir_s)
    }
  }
  
  ## Save the tiles
  for(i in 1:length(list_tiles_names)){
    infile_reg_outline <- list_tiles_names[i]
    tile_spdf <- list_tiles[[i]]
    writeOGR(tile_spdf,dsn= out_dir_s,layer= sub(".shp","",infile_reg_outline), 
             driver="ESRI Shapefile",overwrite_layer="TRUE")
  }
  
  #browser()
  
  list_coord_val <- lapply(list_centroids,function(x)(cbind(coordinates(x)[1],coordinates(x)[2])))
  coord_val <- do.call(rbind,list_coord_val)
  
  ### generate data.frame with ID, filenames and dirname
  df_tiles <- data.frame(ID=1:n_tiles,
                         x_center=coord_val[,1],
                         y_center=coord_val[,2],
                         dirname=out_dir_s,
                         filename=list_tiles_names)
  
  out_filename <- paste0("df_tiles","_",out_suffix,".txt")
  write.table(df_tiles, 
              file=file.path(out_dir_s,out_filename),
              sep=",")
  
  ##### Prepare return object #####
  
  tile_obj <- list(list_tiles,df_tiles)
  names(tile_obj) <- c("list_tiles","df_tiles")
  #View(df_tiles)
  
  return(tile_obj)
}

centroids_shp_fun <- function(i,list_shp_reg_files){
  
  #
  shp_filename <- list_shp_reg_files[[i]]
  layer_name <- sub(".shp","",basename(shp_filename))
  path_to_shp <- dirname(shp_filename)
  shp1 <- try(readOGR(path_to_shp, layer_name)) #use try to resolve error below
  #shp_61.0_-160.0
  #Geographical CRS given to non-conformant data: -186.331747678
  
  #shp1<-readOGR(dirname(list_shp_reg_files[[i]]),sub(".shp","",basename(list_shp_reg_files[[i]])))
  if (!inherits(shp1,"try-error")) {
    pt <- gCentroid(shp1)
    #centroids_pts[[i]] <- pt
  }else{
    pt <- shp1
    #centroids_pts[[i]] <- pt
  }
  
  #shps_tiles[[i]] <- shp1
  #centroids_pts[[i]] <- centroids
  
  shp_obj <- list(shp1,pt)
  names(shp_obj) <- c("spdf","centroid")
  return(shp_obj)
}


########################### end of script #########################
