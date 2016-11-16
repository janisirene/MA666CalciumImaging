function [finalBinaryImage,Components,Centroids] = CaImGetROIs(filename,estNeuronSize,maxNeurons)
%CaImGetROIs.m
%   Take as input a .avi file and get out the individual ROIs that a simple
%    algorithm has identified. Use with a .avi video that Janis created, or
%    it is easy to alter the code to input a .tif file.
%   1) Performs an adaptive Wiener filter on each frame of the image to
%       reduce noise
%   2) Calculates the maximum autocorrelation across time of each pixel
%   3) Performs a statistical test on these maximum autocorrelations to
%       find pixels with significant autocorrelations (accounts for
%       multiple comparisons) to create a binary image
%   4) Performs morphological closing to fill in any gaps for each of the
%       individual ROIs 
%   5) Uses the MATLAB command bwconncomp to identify individual objects in
%       the binary image, get the set of pixels comprising that object, and
%       find object centroids
%
% INPUT: filename - .avi filename as a string
%        estNeuronSize - size of ROI in pixels (defaults to 3)
% OUTPUT: finalBinaryImage - binary image of identified ROIs
%         Components - a structure giving information about the identified
%          ROIs (how many, what pixels do they comprise)
%         Centroids - a structure with each ROI's centroid
%         several figures
% 
%Created: 2016/11/08
% Byron Price
%Updated: 2016/11/10
%  By: Byron Price

% Future steps:
%  1) think about the noise inherent in actual calcium imaging videos
%      e.g. if the noise is correlated in time across more than 1 lag, then
%      many of the noisy pixels will be picked out as having significant
%      autocorrelations
%  2) consider an autocorrelation model that accounts for typical dynamics
%      of calcium fluorescence
%  3) consider what to do if two ROIs overlap (cross-covariance or
%      cross-correlation amongst some window of neighboring pixels)
%  4) test on real data

if nargin < 2
    estNeuronSize = 5;
    maxNeurons = 100;
elseif nargin < 3
    maxNeurons = 100;
end

estNeuronArea = pi*estNeuronSize*estNeuronSize;

if ischar(filename) % data input is .avi
    [~, ~, ext] = fileparts(filename);
    switch ext
        case '.avi'
            vidObj = VideoReader(filename);
            numFrames = vidObj.FrameRate*vidObj.Duration;
            fullVideo = zeros(vidObj.Height,vidObj.Width,numFrames);
            
            count = 1;
            while hasFrame(vidObj)
                temp = readFrame(vidObj);
                fullVideo(:,:,count) = mean(double(temp),3);
                count = count+1;
            end
            fullVideo = fullVideo(60:470,230:640,:);
        case '.tif'
            fullVideo = readTifStack(filename);
            numFrames = size(fullVideo, 3);
    end
else % data input is already 3d array
    fullVideo = filename;
    numFrames = size(fullVideo, 3);
end

% denoise with wiener filter
width = size(fullVideo,1);
height = size(fullVideo,2);

fltVideo = zeros(size(fullVideo));
for ii=1:numFrames
    temp = fullVideo(:,:,ii);
    fltVideo(:,:,ii) = wiener2(temp,[5,5]);
end

%implay(uint8(fullVideo));


% obtain autocorrelation image
autoCorrImg = zeros(width,height);
for ii=1:width
    for jj=1:height
        [acf,~] = autocorr(squeeze(fltVideo(ii,jj,:)));
        % pixels with a temporal structure will have larger acf at nonzero
        % lags
        if isnan(acf) ~= 1
            autoCorrImg(ii,jj) = max(abs(acf(acf<1)));
        end
    end
end

figure();imagesc(autoCorrImg);caxis([0 1]);colormap('jet');colorbar;
title('Maximum Autocorrelation Image');


% for thresholding ... perform statistical test on autocorrelation
%  coefficients
numComparisons = width*height;
% alpha = 1-0.05/numComparisons;
% 
% q = norminv(alpha,0,1);
% threshold = q/sqrt(numFrames); 

totalNeuronArea = maxNeurons*estNeuronArea;
threshold = quantile(autoCorrImg(:),1-totalNeuronArea/numComparisons);
binaryAutoCorr = autoCorrImg > threshold;

figure();histogram(autoCorrImg(:));hold on;
[N,~] = histcounts(autoCorrImg(:));
plot(ones(100,1).*threshold,linspace(0,max(N),100),'r','LineWidth',2);
xlabel('Maximum Autocorrelation Coefficient');ylabel('Count');
title('Histogram of Maximum Autocorrelation Coefficients');
legend('Histogram','Bonferroni-Corrected Threshold');

% morphological opening and closing?
se = strel('disk',round(estNeuronSize/4));
se2 = strel('disk',round(estNeuronSize));
finalBinaryImage = imclose(imopen(binaryAutoCorr,se),se2);

figure();imagesc(finalBinaryImage);colormap('bone');
title('Binary Mask for ROI Detection');

% remove noisy background
maskedVideo = zeros(size(fltVideo));
for ii=1:numFrames 
    maskedVideo(:,:,ii) = fltVideo(:,:,ii).*finalBinaryImage;
end

maxlag = 5;
%% summed crosscorr
%{
% look for components to either merge or separate with cross-correlation
summedCrossCorr = zeros(width,height);
divisor = zeros(width,height);
for ii=1:width
    for jj=1:height
            for kk=-2:2
                for ll=-2:2
                    if (ii+kk) > 0 && (jj+ll) > 0 && (ii+kk) <= width && (jj+ll) <= height %&& kk ~= 0 && ll ~= 0
                        summedCrossCorr(ii+kk,jj+ll) = summedCrossCorr(ii+kk,jj+ll)+max(xcorr(squeeze(fltVideo(ii,jj,:)),squeeze(fltVideo(ii+kk,jj+ll,:)),maxlag,'coeff'));
                        divisor(ii+kk,jj+ll) = divisor(ii+kk,jj+ll)+1;
                    end
                end
            end
    end
end
summedCrossCorr = summedCrossCorr./divisor;
figure();imagesc(summedCrossCorr);title('Summed Cross-Correlation Image');

% at this point, we could try either a non-parametric statistical test or
%  attempt to figure out the distribution of these coefficients and then 
%  eliminate regions in the image with low summed cross-correlation
%  coefficients (indicating that they lie at the border between two cells)
summedCrossCorr = summedCrossCorr(:);
figure();histogram(summedCrossCorr(summedCrossCorr~=0));
title('Histogram of Summed Cross-Correlation Coefficients');
xlabel('Coefficient Magnitude');ylabel('Count');
%}

%% adjacency matrix and hierarchical clustering
% store values in a sparse matrix by making an index array and a value
% array
usePixels = find(finalBinaryImage); % only use pixels deemed signal
[tempr, tempc] = meshgrid(1:length(usePixels), 1:length(usePixels));
indexArray = [tempr(:), tempc(:)]; % index of usePixel
kill = (indexArray(:, 1) <= indexArray(:, 2));  % unique pairs
indexArray(kill, :) = [];
usePixelArray = usePixels(indexArray); % actual pixel indices in image
clear tempr tempc;

% distance between pixels - skip pairs that are too far apart or pairs that
% are the same pixel (or leave them all in to use linkage)
col = ceil(usePixelArray / width);
row = usePixelArray - width * (col - 1);
pixelDist = sqrt((col(:, 1) - col(:, 2)).^2 + (row(:, 1) - row(:, 2)).^2);

xcorrArray = zeros(length(indexArray), 1);
for ii = 1:length(xcorrArray)
    if pixelDist(ii) > 2 * estNeuronSize
        continue; 
    end
    idxr = row(ii, 1);
    idxc = col(ii, 1);
    jdxr = row(ii, 2);
    jdxc = row(ii, 2);

    xcorrArray(ii) = max(xcorr(squeeze(maskedVideo(idxr, idxc, :)),...
        squeeze(maskedVideo(jdxr, jdxc, :)), maxlag, 'unbiased'));
end
xcorrArray = max(0, xcorrArray); % less than 0 just ignore it
metric = max(xcorrArray) - xcorrArray; % 
Z = linkage(metric', 'average');
t = cluster(Z);
keyboard;


% bwconncomp finds groups in the binary image
Components = bwconncomp(finalBinaryImage);
Centroids = regionprops(Components,'Centroid');
end

