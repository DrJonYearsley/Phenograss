# Script to split MODIS data up into big blocks for processing
#
# Jon Yearsley
# Jon.Yearsley@ucd.ie
# June 2022
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rm(list=ls())

dataDir = '/media/jon/MODIS_data/'
outputDir = '~/Research/Data/MODIS_pasture/'
yearStr = paste0('A',c(2003:2022))

# Load coordinates
load(file.path(dataDir,'modis_pasture_coords.RData'))



x_values=round(unique(coords_df$x),1)
y_values=round(unique(coords_df$y),1)

x_blocks = seq(from=min(x_values), to=max(x_values)+1, length.out=4)
y_blocks = seq(from=min(y_values), to=max(y_values)+1, length.out=6)


for (y in yearStr) {
  filein = file.path(dataDir,'MODIS_AllIreland',paste0('modis_pasture_',y,'.RData'))
  load(filein)
  
  print(filein)

  block=1
  for (i in 2:length(x_blocks)) {
    ind_x = d$x_MODIS>=x_blocks[i-1] & d$x_MODIS<x_blocks[i]
    for (j in 2:length(y_blocks)) {
      ind_y = d$y_MODIS>=y_blocks[j-1] & d$y_MODIS<y_blocks[j]
      
      if (any(ind_x & ind_y)) {
        d_block = d[ind_x & ind_y,]
        d_block$pixelID = paste0('x',match(round(d_block$x_MODIS,1),x_values),
                                 'y',match(round(d_block$y_MODIS,1),y_values))
        
        if (block==1) {
          print(head(d_block))
        }
        
        fileOut = file.path(outputDir,paste0('MODIS_block_',block,'_',y,'.RData'))
        save('d_block',file=fileOut)
      }
      block = block+1
    }
  }
}