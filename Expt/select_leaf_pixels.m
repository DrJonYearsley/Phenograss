function rgb_leaves = select_leaf_pixels(im)
% Function takes a demosaiced and color corrected image 
% and extracts pixels that correspond to leaves in the image
%
% Inputs:
%    im  - input image
%
% Outputs:
%    rgb_leaves  - array of RGB values for each pixel identified as being
%    part of a leaf
%
% Jon Yearsley (jon.yearsley@ucd.ie)
% Nov 2021
% ++++++++++++++++++++++++++++++++++++++++++++++++++++


%% Convert image to HSV color space
% 
% Then identify the paper on which the leaves are placed

I = rgb2hsv(im);
sz = size(I);

% Identify paper as having a value>=0.75
mask_paper = I(:,:,3) >= 0.75;
[xInd, yInd] = find(mask_paper==1);
[a,b] = hist(xInd,sz(1));
xLims = round([min(b(a>max(a)/2)), max(b(a>max(a)/2))]);
[a,b] = hist(yInd,sz(2));
yLims = round([min(b(a>max(a)/2)), max(b(a>max(a)/2))]);
im_paper = im(xLims(1):xLims(2), yLims(1):yLims(2),:);
I_paper = I(xLims(1):xLims(2), yLims(1):yLims(2),:);

clear("I","im","mask_paper")

%% Select colors that are in the yellow/green/blue part of the spectrum

% A range of hues to select
channel1Min1 = 0.000;
channel1Max1 = 0.400;

channel1Min2 = 0.875;
channel1Max2 = 1.000;

% A range of saturations
channel2Min = 0.200;
channel2Max = 1.000;

% A range of values
channel3Min = 0.000;
channel3Max = 0.900;

% Apply the mask
mask1a = (I_paper(:,:,1) >= channel1Min1 ) & (I_paper(:,:,1) <= channel1Max1) & ...
   (I_paper(:,:,2) >= channel2Min ) & (I_paper(:,:,2) <= channel2Max) & ...
   (I_paper(:,:,3) >= channel3Min ) & (I_paper(:,:,3) <= channel3Max);

mask1b = (I_paper(:,:,1) >= channel1Min2 ) & (I_paper(:,:,1) <= channel1Max2) & ...
   (I_paper(:,:,2) >= channel2Min ) & (I_paper(:,:,2) <= channel2Max) & ...
   (I_paper(:,:,3) >= channel3Min ) & (I_paper(:,:,3) <= channel3Max);

mask1 = mask1a | mask1b;

% Perform some adaptive pixel selection starting from mask1
mask2 = bwareaopen(mask1,100);
mask3 = activecontour(im_paper,mask2,50);

% Finally impose another color filter (filters out low saturation)
channel1MinA = 0.000;
channel1MaxA = 1.000;

channel2MinA = 0.250;
channel2MaxA = 1.000;

channel3MinA = 0.000;
channel3MaxA = 1.000;

mask4 = (I_paper(:,:,1) >= channel1MinA ) & (I_paper(:,:,1) <= channel1MaxA) & ...
   (I_paper(:,:,2) >= channel2MinA ) & (I_paper(:,:,2) <= channel2MaxA) & ...
   (I_paper(:,:,3) >= channel3MinA ) & (I_paper(:,:,3) <= channel3MaxA);

mask_final = mask4 & mask3;


im_final = im_paper.*repmat(uint16(mask_final),1,1,3);

clear("mask1a","mask1b","mask1","mask2","mask3","mask4")

%% Extract RGB values from the identified pixels

% Count pixels in mask
nPixel = sum(mask_final(:)==1);

% Array to store RGB values
rgb_leaves = zeros(nPixel,3);

% Loop through the three channels (RGB) and extract the values
for i=1:3
    tmp = im_final(:,:,i);
    rgb_leaves(:,i) = tmp(mask_final);
end

