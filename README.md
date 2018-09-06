Generating Remote processing workflow with specific example of application to GIMMS.

The implementation of the workflow is aimed at being general. Scripts are still under
development but here are the main steps with associated scripts for documentation:

1) Get study area
https://github.com/bparment1/spatial_variation_GIMMS/blob/master/get_study_region_data.R

This set of functions generates an area based on a country/region using GADM (https://gadm.org/) openly available dataset.

2) Split area in tiles

https://github.com/bparment1/spatial_variation_GIMMS/blob/master/generate_spatial_tiles_functions.R

This set of functions generates tiles with various options including number of tiles, overlapping % etc. The gnerated tile system is used in processing the study area tile by tile.

3) Specific apply function for time series in a tile

https://github.com/bparment1/spatial_variation_GIMMS/blob/master/lag_processing_functions.R

In this exampe, I use functions to extract Moran's I lag of images in time series. The function can be substituted for other analyses,

4) Mosaic back to obtain modeled results for a given area

Tiles with model results are mosaiced back together. This part can be expensive in term of computation so 
I am modifying earlier code using gdal_merge.py modified for specific use (by Alberto Guzman). 
In case of overlapping tiles, the tiles must be weighted so that there are no edge effect. 
The current code implements various options for the average weighting managing the large number of files produce during processing.

