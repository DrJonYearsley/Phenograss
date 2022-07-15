% PROCESS_HARVEST_IMAGES.M
%
% Process ryegrass leaf images from Phenograss experiment
%
% Each image has a colour swatch that is used to correct the colours in the
% image.
%
% Run this script after navigating in MATLAB to the folder that holds the images
% Matlab is expecting these images to be RAW images with a suffix .RW2
% The script will process all images in the folder with a .RW2 suffix
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
        outputFile = join([outputFileBase{1} "tiff"],".");
        imwrite(im, outputFile,'tiff','Compression','none','ColorSpace','rgb')
    end
    % Extract pixels from leaves
    rgb_leaves = select_leaf_pixels(im);
    
    leafOutputFile = join([outputFileBase{1} "leaf_RGB_pixels.mat"],"_");
    save(leafOutputFile, 'rgb_leaves')
end


