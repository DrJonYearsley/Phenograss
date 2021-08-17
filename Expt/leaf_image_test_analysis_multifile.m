% LEAF_IMAGE_TEST_ANALYSIS_MULTIFILE.M
%
% Process ryegrass leaf images from Phenograss experiment
%
% Each image has a colour swatch that is used to correct the colours in the
% image.
%
% The script is the same as LEAF_IMAGE_TEST_ANALYSIS.m but it
% iterates over multiple files
%
% Steps:
%   0. List files in directory which are to be processed
%   1. Import image
%   2. Identify colour swatch
%   3. Correct colours in whole image
%   4. Identify white paper
%   5. Identify leaves in region of white paper (based on HSV thresholding)
%   6. Extract RGB values
%   7. Display a histogram of R, G and B channels
%   8. Repeat steps 1-7 for each file
%
% Jon Yearsley (Jon.Yearsley@ucd.ie
% Aug 2021
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++


clear all

display = false;   % If true display graphs and images

%% Find files to import
files = dir('*.jpg');
nFiles = numel(files);



%% Loop through all files
for f=1:nFiles
    %% Import image and correct colours
    inputFile = files(f).name;
    RGB = imread(inputFile);
    
    disp(['Processing file ', inputFile])
    
    
    % Identify the color cahrt within the image
    chart = colorChecker(RGB);
    
    
    % Compare measured colour to reference colour values and calculate colour
    % correction matrix (ccm)
    [colorTable, ccm] = measureColor(chart);
    
    
    % Correct colours in the image
    RGB_cc = imapplymatrix(ccm(1:3,:)',RGB,ccm(4,:));

    
    if display
        % Display original image
        imshow(RGB)
        
        % Display the colour chart
        displayChart(chart)
        
        % Display a comparison of reference and measured colours
        displayColorPatch(colorTable)
        
        
        % Display corrected image
        imshow(RGB_cc)
    end
    
    
    %% Convert RGB image to HSV color space
    I = rgb2hsv(RGB_cc);
    
    sz = size(I);
    
    
    %% Find the paper limits
    
    % Threshold based on values greater than 0.75
    mask_paper = I(:,:,3) >= 0.75;
    
    [xInd, yInd] = find(mask_paper==1);
    
    % Then look for rows and columns with lots of values>0.75 (i.e. the sheet
    % of paper)
    [a,b] = hist(xInd,sz(1));
    xLims = round([min(b(a>max(a)/2)), max(b(a>max(a)/2))]);
    
    [a,b] = hist(yInd,sz(2));
    yLims = round([min(b(a>max(a)/2)), max(b(a>max(a)/2))]);
    
    % Extract out the paper part of the image
    RGB_paper = RGB_cc(xLims(1):xLims(2), yLims(1):yLims(2),:);
    I_paper = I(xLims(1):xLims(2), yLims(1):yLims(2),:);
    
    if display
        % Display the cropped image
        imshow(RGB_paper)
    end
    
    %% Segment the paper part of the image
    
    
    % Define hue mask (select out green tones)
    channel1Min = 0.144;
    channel1Max = 0.514;
    
    
    % Define saturation mask
    %(remove low saturation values, which are close to neutral)
    channel2Min = 0.1;
    channel2Max = 1.000;
    
    % Define thresholds for value setting
    % (low values are dark, high values are light)
    channel3Min = 0.000;
    channel3Max = 0.5;
    
    % Create mask based on chosen histogram thresholds
    mask = (I_paper(:,:,1) >= channel1Min ) & (I_paper(:,:,1) <= channel1Max) & ...
        (I_paper(:,:,2) >= channel2Min ) & (I_paper(:,:,2) <= channel2Max) & ...
        (I_paper(:,:,3) >= channel3Min ) & (I_paper(:,:,3) <= channel3Max);
    
    if display
        % visualise the mask
        imshow(mask)
        imshow(RGB_paper.*repmat(uint8(mask),1,1,3))
    end
    
    %% Obtain red, green and blue values after mask
    
    % Count pixels in mask
    nPixel = sum(mask(:)==1);
    image_sz = size(RGB_paper);
    
    % Find indices of mask (used to keep spatial information about pixels)
    [row_pixelInd, col_pixelInd] = find(mask);
 
    
    % Array to store RGB values
    rgb_leaves = zeros(nPixel,3);
    
    labs=["Red values", "Green values","Blue Values"];
    for i=1:3
        ind = sub2ind(image_sz, row_pixelInd, col_pixelInd, i*ones(nPixel,1));
        rgb_leaves(:,i) = RGB_paper(ind);
        
        if display
            % Display histograms of rgb values from leaves
            subplot(3,1,i)
            hist(rgb_leaves(:,i))
            xlabel(labs(i))
        end
    end
    
    %% Save data for this file
    filename = ['processed_leaf_data_file', num2str(f), '.mat'];
    save(filename,'rgb_leaves','nPixel','inputFile','row_pixelInd','col_pixelInd')
    
end
%% Perform PCA of colour values

[coeff,score,latent,tsquared,explained,mu] = pca(rgb_leaves);

