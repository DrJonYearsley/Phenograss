# Processing of Digital Images from the Growth Chamber

This folder contains 3 Matlab scripts that perform 
The image processing from RAW image files. 



## process_harvest_images.m
      This is the main Matlab script that calls the other scripts.



## process_raw_image.m
	This script performs the demosaicing of the RAW image and the color correction based upon the standard coloraturas swatch in the image



## select_leaf_pixels.m
	This script takes the color corrected image and creates a mask that selects only pixels that are identified as being a leaf