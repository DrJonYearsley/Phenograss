function srgb_image2 = process_raw_image(filename, method)
    arguments
        filename char
        method   double = 1
    end
% Function imports a RAW image file, corrects the image, demosaics
% and finally uses a color chart to perform a final correction. 
%
% The pipeline is based upon  
%  https://uk.mathworks.com/help/images/end-to-end-implementation-of-digital-camera-processing-pipeline.html
%
% Inputs:
%     filename   -   the RAW image to be processed
%     method     -   color correction method
%                    1 = (default)
%                    2 = use color matrix from the camera
%
% Outputs:
%     srgb_image1  - processed RGB image 
%
%
% Jon Yearsley (jon.yearsley@cd.ie)
% Nov 2021
% =======================================


%% Read the raw image (a raw colour filter array, CFA)
cfaImage = rawread(filename);
cfaInfo = rawinfo(filename);


%% Display some info
colorInfo = cfaInfo.ColorInfo;

%% Scale pixel values to a suitable range

% Perform Black Level Correction
% RAW images do not have a true black value. Even with the shutter closed, 
% electricity flowing through the sensors causes nonzero photon counts. 
% Cameras use the value of the masked pixels to compute the black level of 
% the CFA image. To scale the image, subtract the measured black level from the CFA image data.

% RAW file formats can report this black level in different formats. 
% The RAW file metadata for the image in this example specifies black level 
% as a vector, with one element per channel of the CFA image. Other RAW 
% file formats, such as DNG, specify black level as a repeated m-by-n matrix 
% that starts at the top left corner of the visible portion of the CFA.


% Get the black level value of the RAW data from the BlackLevel metadata field.
blackLevel = colorInfo.BlackLevel;

% To perform black level correction, first convert the black level vector to a 2-by-2 matrix. 
% Omit this step for RAW images that specify black level as a matrix.
blackLevel = reshape(blackLevel,[1 1 numel(blackLevel)]);
blackLevel = planar2raw(blackLevel);

% Replicate the black level matrix to be the size of the visible image.
repeatDims = cfaInfo.ImageSizeInfo.VisibleImageSize ./ size(blackLevel);
blackLevel = repmat(blackLevel,repeatDims);

% Subtract the black level matrix from the CFA image matrix.
cfaImage = cfaImage - blackLevel;

% ==========================
% Clamp negative pixel values

% To correct for CFA data values less than the black-level value, clamp the values to 0.
cfaImage = max(0,cfaImage);


% ==========================
% Scale pixel values

% RAW file metadata often represent the white level as the maximum value 
% allowed by the data type. If this white level value is much higher than 
% the highest intensity value in the image, then using this white level 
% value for scaling results in an image that is darker than it should be. 
% To avoid this, scale the CFA image using the maximum pixel value found in the image.
cfaImage = double(cfaImage);
maxValue = max(cfaImage(:));

cfaImage = cfaImage ./ maxValue;


%% Adjust White Balance
% White balance is the process of removing unrealistic color casts from a 
% rendered image, such that it appears closer to how human eyes would see the subject.

% Get the white balance values from the metadata. There are two types of 
% white balance metadata available. This step of the example uses the 
% CameraAsTakenWhiteBalance field scales the color channels to balance 
% the linear pixel values. The example uses the D65WhiteBalance field later 
% in the pipeline to adjust the colors to a D65 white point.

whiteBalance = colorInfo.CameraAsTakenWhiteBalance;

% Scale the multipliers so that the values of the green color channels are 1.
gLoc = strfind(cfaInfo.CFALayout,"G"); 
gLoc = gLoc(1);
whiteBalance = whiteBalance/whiteBalance(gLoc);

whiteBalance = reshape(whiteBalance,[1 1 numel(whiteBalance)]);
whiteBalance = planar2raw(whiteBalance);

% Replicate the white balance matrix to be the size of the visible image.
whiteBalance = repmat(whiteBalance,repeatDims);
cfaWB = cfaImage .* whiteBalance;

% Convert the CFA image to a 16-bit image.
cfaWB = im2uint16(cfaWB);

%% Demosaic
% Convert the Bayer-encoded CFA image into a truecolor image by demosaicing. 
% The truecolor image is in linear camera space.

cfaLayout = cfaInfo.CFALayout;
imDebayered = demosaic(cfaWB,cfaLayout);

%% Two methods to convert to RGB space

if method==1
    % Method 1: Convert from Camera Color Space to RGB Color Space

    % Get the transformation matrix between the linear camera space and the
    % XYZ profile connection space from the CameraToXYZ metadata field. 
    % This matrix imposes an RGB order.
    cam2xyzMat = colorInfo.CameraToXYZ;

    % Normalize the cam2xyzMat matrix according to a D65 white point.
    % Get the XYZ normalization values from the D65WhiteBalance metadata field.
    whiteBalanceD65 = colorInfo.D65WhiteBalance;

    % The white balance multipliers are ordered according to the CFALayout metadata field.
    % Reorder the multipliers to match the row ordering of the cam2xyzMat matrix.
    cfaLayout = cfaInfo.CFALayout;
    wbIdx(1) = strfind(cfaLayout,"R");
    gidx = strfind(cfaLayout,"G");
    wbIdx(2) = gidx(1);
    wbIdx(3) = strfind(cfaLayout,"B");

    wbCoeffs = whiteBalanceD65(wbIdx);
    cam2xyzMat = cam2xyzMat ./ wbCoeffs;

    % Convert the image from the linear camera space to the XYZ color space
    % using the imapplymatrix function. Then, convert the image to the sRGB
    % color space and apply gamma correction using the xyz2rgb function.
    imXYZ = imapplymatrix(cam2xyzMat,im2double(imDebayered));
    srgb_image1 = xyz2rgb(imXYZ,"OutputType","uint16");


else
    % Method 2: Convert from Camera Color Space to RGB Color Space

    % Use Conversion Matrix from RAW File Metadata
    % Convert the image from the linear camera space to the linear RGB color
    % space using the transformation matrix in the CameraTosRGB metadata field.
    cam2srgbMat = colorInfo.CameraTosRGB;
    imTransform = imapplymatrix(cam2srgbMat,imDebayered,"uint16");

    % Apply gamma correction to bring the image from the linear sRGB color space
    % to the sRGB color space.
    srgb_image1 = lin2rgb(imTransform);
end


%% Final color correction using the X-rite colorchart classic

% Method 1
chart = colorChecker(srgb_image1, sensitivity=0.8);
[~, ccm] = measureColor(chart);
srgb_image2 = imapplymatrix(ccm(1:3,:)',srgb_image1,ccm(4,:));
