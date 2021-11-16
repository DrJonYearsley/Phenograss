% PROCESS_HARVEST_IMAGES.M
%
% Process ryegrass leaf images from Phenograss experiment
%
% Each image has a colour swatch that is used to correct the colours in the
% image.
%
%
% Steps:
%   0. List files in directory which are to be processed
%   1. Use process_raw_image() to demosaic RAW image
%   2. Save demosaiced image as a tiff
%   3. Use select_leaf_pixels() to color correct and identify leaf pixels in the image
%       a. Identify colour swatch
%       b. Correct colours in whole image
%       c. Identify white paper
%       d. Identify leaves in region of white paper (based on HSV thresholding)
%       e. Extract RGB values
%    4. Save leaf pixel data
%
% Jon Yearsley (Jon.Yearsley@ucd.ie)
% Aug 2021
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++


clear all
clear functions

% If true then save demosaiced image as a tiff
save_tiff = true;


%% Find files to import
filenames = dir('*.RW2');
nFiles = numel(filenames);



%% Extract and color correct image

for f=1:nFiles
    % Import raw image and demosaic
    inputFile = filenames(f).name;
    im = process_raw_image(inputFile, 2);

    % Write image as a tiff (with no compression and an RGB colorspace)
    outputFileBase = split(inputFile,".RW2");
    if (save_tiff) 
      imwrite(im, outputFileBase{1},'tiff','Compression','none','ColorSpace','rgb')
    end
    % Extract pixels from leaves
    rgb_leaves = select_leaf_pixels(im);
    
    leafOutputFile = join([outputFileBase{1} "leaf_RGB_pixels.mat"],"_");
    save(leafOutputFile, 'rgb_leaves')
end



% %% Visualise the data
% figure(1)
% imshow(im2)
% 
% %%
% figure(2)
% 
% max_val = 2^16-1;
% width=500;
% 
% subplot(3,1,1)
% hist(rgb_leaves(:,1), round(max_val/width))
% xlim([0, 6e4])
% subplot(3,1,2)
% hist(rgb_leaves(:,2), round(max_val/width))
% xlim([0, 6e4])
% subplot(3,1,3)
% hist(rgb_leaves(:,3), round(max_val/width))
% xlim([0, 6e4])
