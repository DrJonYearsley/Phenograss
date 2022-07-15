# Remote sensing analysis of MODIS data

This folder contains the R scripts that perform 
the analysis of the remote-sensing data



**create_quadrat_polygons.R**
This creates shape files for the 10 km square quadrats used in the analysis


**land_cover_MODIS.R**
This creates a raster giving the % of pasture in each MODIS pixel across Ireland



**phenograss_process_MODIS_v2.R**
This ingests the raw MODIS data and produces cleaned data for each 10 km square


**segmentation_parrallel.R**
This performs the segmentation analysis which estimates the phenophases

**segmentation_functions.R**
Contains the basic functions for the segmentation analysis


**analyse_MODIS_phenology_finalreport.R**
This performs a spatial analysis of the estimated phenophases
