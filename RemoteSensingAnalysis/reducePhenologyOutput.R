# Simplify the phenology data for analysis
#
#
# Jon Yearsley (jon.yearsley@ucd.ie)
# March 2022
# ++++++++++++++++++++++++++++++++++++++++++++++

rm(list=ls())



dataDir = '/media/jon/MODIS_data/PhenologyOutput_test/'
outputDir = '/media/jon/MODIS_data/PhenologyOutput_toAnalyse/'

input_file_preffix = 'phenology'

# Import data --------
squares = c(1:21)

colsToKeep = c("pixelID","square","x_ITM","y_ITM","x_MODIS","y_MODIS",
               "year","phase","t","lowerCI","upperCI","warning")

# Filename segmented data
for (i in 1:length(squares)) {
  print(paste("Importing data for square",squares[i],":",Sys.time()))
  
  files = list.files(path=dataDir, 
                     pattern=paste0(input_file_preffix,"_square_",squares[i],"_"),
                     full.names=TRUE)
  
  for (f in 1:length(files)) {
    load(files[f])
    
    
    ind = segment_output$model=="smooth" & 
      !segment_output$warning & 
      segment_output$phase%in%c(1,2,3)
    if (f==1 ) {
      phenology = segment_output[ind,colsToKeep]
    } else {
      phenology = rbind(phenology,
                        segment_output[ind,colsToKeep])
    }
    rm(list=c("segment_output","d_final","knots","year","nSegBreaks","s"))
  }
  
  filename = paste0("phenology_square_",squares[i],'.RData')
  save(phenology, file = file.path(outputDir, filename))
  
  rm(list=c("phenology"))
  
}
