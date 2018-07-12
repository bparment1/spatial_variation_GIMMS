
  

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
  
  selected_countries_spdf <-getData("GADM", 
                               country=codes_countries$ISO3[index_values]
                               ,level=0)
  
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


country_names <- c("Kazakhstan") #can have mulitple names as a vector here
#undebug(get_countries_outlines)
reg_outline_spdf <- get_countries_outlines(country_names,
                       out_dir=".",
                       out_suffix="kazakhstan")
#
plot(reg_outline_spdf)  

country_names <- c("Puerto Rico") #can have mulitple names as a vector here
#undebug(get_countries_outlines)
reg_outline_spdf <- get_countries_outlines(country_names,
                                           out_dir=".",
                                           out_suffix="puerto_rico")


#Select tiles matching the areas of interest and with specific content:
#Kenya<-getData("GADM", country="KE", level=0)
KZ_spdf <-getData("GADM", country="KZ", level=0)
KZ_sf <- as(KZ_spdf,"sf")

plot(KZ_sf$geometry,add=T,border="blue")
#Kenya1<-getData("GADM", country="KE", level=1)

st_crs(KZ_sf)

selected_poly_ID <- st_intersects(KZ_sf,tiles_combined_sf)
selected_poly_ID <- unlist(selected_poly_ID)

tiles_reg_sf <- filter(tiles_combined_sf,ID%in%selected_poly_ID)

### visualize selected polygons
plot(r_ref)
text(df_tiles$x_center,df_tiles$y_center,df_tiles$ID,cex=0.5)
plot(KZ_sf$geometry,add=T,border="blue")
plot(tiles_reg_sf$geometry,border="red",add=T)

## Let's drop 41 since the overlap is low

tiles_reg_sf <- filter(tiles_reg_sf,ID!=41)

out_tiles_filename <- file.path(out_dir,"tiles_reg.shp")
st_write(tiles_reg_sf,
         dsn=out_tiles_filename)
