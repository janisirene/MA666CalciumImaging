function [finalBinaryImage,Components,Centroids] = CaImGetROIs(filename,estNeuronSize,maxNeurons)
%CaImGetROIs.m
%   Take as input a .avi or .tif file and get out the individual ROIs that a simple
%    algorithm has identified.
%   1) Performs an adaptive Wiener filter on each frame of the image to
%       reduce noise
%   2) Calculates a spatial average maximum cross-correlation across time 
%       at each pixel (maximum crosss-correlation with small collection of
%       neighboring pixels)
%   3) Performs a statistical test on these cross-correlations to
%       find pixels with significant spatiotemporal correlations with itself
%       and it's neighbors to create a binary image
%   4) Performs morphological closing to fill in any gaps for each of the
%       individual ROIs and opening to remove very small ROIs
%   5) Clustering to merge or separate tentative ROIs
%   5) Uses the MATLAB command bwconncomp to identify individual objects in
%       the binary image, get the set of pixels comprising that object, and
%       find object centroids
%
% INPUT: filename - .avi or .tif filename as a string
%        estNeuronSize - size of ROI in pixels (defaults to 3)
%        maxNeurons - maximum number of neurons in field of view (defaults
%         to 100) 
% OUTPUT: finalBinaryImage - binary image of identified ROIs
%         Components - a structure giving information about the identified
%          ROIs (how many, what pixels do they comprise)
%         Centroids - a structure with each ROI's centroid
%         several figures
% 
%Created: 2016/11/08
% Byron Price
%Updated: 2016/11/15 
%  By: Byron Price & Janis Intoy

% Future steps:
%  1) think about the noise inherent in actual calcium imaging videos
%      e.g. if the noise is correlated in time across more than 1 lag, then
%      many of the noisy pixels will be picked out as having significant
%      autocorrelations
%  2) clustering algorithms
%  3) remove non-stationary background 

if nargin < 2
    estNeuronSize = 5;
    maxNeurons = 100;
elseif nargin < 3
    maxNeurons = 100;
end

estNeuronSize = round(estNeuronSize);
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
    fltVideo(:,:,ii) = wiener2(temp,[estNeuronSize,estNeuronSize]);
end

%implay(uint8(fullVideo));


% % obtain autocorrelation image
% autoCorrImg = zeros(width,height);
% for ii=1:width
%     for jj=1:height
%         [acf,~] = autocorr(squeeze(fltVideo(ii,jj,:)));
%         % pixels with a temporal structure will have larger acf at nonzero
%         % lags
%         if isnan(acf) ~= 1
%             autoCorrImg(ii,jj) = max(abs(acf(acf<1)));
%         end
%     end
% end
% 
% figure();imagesc(autoCorrImg);caxis([0 1]);colormap('jet');colorbar;
% title('Maximum Autocorrelation Image');

% look for components to either merge or separate with cross-correlation
summedCrossCorr = zeros(width,height);
divisor = zeros(width,height);
maxlag = 5;
for ii=1:width
    for jj=1:height
            for kk=-estNeuronSize:estNeuronSize
                for ll=-estNeuronSize:estNeuronSize
                    if (ii+kk) > 0 && (jj+ll) > 0 && (ii+kk) <= width && (jj+ll) <= height %&& kk ~= 0 && ll ~= 0
                        summedCrossCorr(ii+kk,jj+ll) = summedCrossCorr(ii+kk,jj+ll)+max(xcorr(squeeze(fltVideo(ii,jj,:)),squeeze(maskedVideo(ii+kk,jj+ll,:)),maxlag,'coeff'));
                        divisor(ii+kk,jj+ll) = divisor(ii+kk,jj+ll)+1;
                    end
                end
            end
    end
end
summedCrossCorr = summedCrossCorr./divisor;
figure();imagesc(summedCrossCorr);title('Spatial Average Maximum Cross-Correlation Image');

% for thresholding ... perform statistical test on autocorrelation
%  coefficients
numComparisons = width*height;
% alpha = 1-0.05/numComparisons;
% 
% q = norminv(alpha,0,1);
% threshold = q/sqrt(numFrames); 

totalNeuronArea = maxNeurons*estNeuronArea;
threshold = quantile(summedCrossCorr(:),1-totalNeuronArea/numComparisons);
binaryCrossCorr = summedCrossCorr > threshold;

forHistogram = binaryCrossCorr(:);
forHistogram = forHistogram(forHistogram ~= 0);
figure();histogram(forHistogram);hold on;
[N,~] = histcounts(forHistogram);
plot(ones(100,1).*threshold,linspace(0,max(N),100),'r','LineWidth',2);
xlabel('Spatial Average Maximum Cross-Correlation');ylabel('Count');
title('Histogram of Maximum Cross-Correlation Coefficients');
legend('Histogram','Bonferroni-Corrected Threshold');



% morphological opening and closing?
se = strel('disk',round(estNeuronSize/4));
se2 = strel('disk',round(estNeuronSize));
finalBinaryImage = imclose(imopen(binaryCrossCorr,se),se2);

figure();imagesc(finalBinaryImage);colormap('bone');
title('Binary Mask for ROI Detection');

% remove noisy background
maskedVideo = zeros(size(fltVideo));
for ii=1:numFrames 
    maskedVideo(:,:,ii) = fltVideo(:,:,ii).*finalBinaryImage;
end

% h = fspecial('laplacian');
% figure();imagesc(filter2(h,summedCrossCorr,'same'));
% at this point, we could try either a non-parametric statistical test or
%  attempt to figure out the distribution of these coefficients and then 
%  eliminate regions in the image with low summed cross-correlation
%  coefficients (indicating that they lie at the border between two cells)


% bwconncomp finds groups in the binary image
Components = bwconncomp(finalBinaryImage);
Centroids = regionprops(Components,'Centroid');
end
