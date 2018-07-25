############### General utility function ########## 
## Help in processing modis data for the workshop group.
##
##
## DATE CREATED: 07/06/2018
## DATE MODIFIED: 07/25/2018
## AUTHORS: Benoit Parmentier  
## Version: 1
## PROJECT: General utility function
## ISSUE: 
## TO DO: - add for other entities (counties, metropolitan areas) etc.
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
  

get_countries_outlines <- function(country_names,out_dir,out_suffix){
  # General function to obtain country outlines in vector format
  #
  #
  ## Add option for geojson?
  
  ####### Start script #####
  
  #test <-getData("GADM",country="*", level=0)
  codes_countries <- ccodes()
  #test <- getData("GDAM",codes_countries$ISO3,level=0)
  codes_countries$ISO3
  ### can be more than one
  #index_values <- which(codes_countries$NAME=="Puerto Rico")
  index_values <- which(codes_countries$NAME%in%country_names)
  
  selected_countries_spdf <- raster::getData("GADM", 
                               country=codes_countries$ISO3[index_values],
                               level=0)
  
  #plot(selected_countries_spdf)
  #class(selected_countries)
  
  outfile_reg_outline <- paste("reg_out_line_",out_suffix,".shp",sep="")
  writeOGR(selected_countries_spdf,
           dsn= out_dir,
           layer= sub(".shp","",outfile_reg_outline), 
           driver="ESRI Shapefile",
           overwrite_layer="TRUE")
  
  return(selected_countries_spdf)
  
}


##################################### End of script ######################