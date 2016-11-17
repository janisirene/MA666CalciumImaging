function [tempBinaryImage,Components,Centroids] = CaImGetROIs(filename,estNeuronRadius,maxNeurons)
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
% INPUT: filename - .avi or .tif filename as a string, or a 3d array
%        estNeuronRadius - radius of ROI in pixels (defaults to 3)
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
%
% TO DO: instead of xcorr do cross correlatin calcs in fourier domain
% (ifft of cross spectrum between time courses of two pixels)

% Future steps:
%  1) think about the noise inherent in actual calcium imaging videos
%      e.g. if the noise is correlated in time across more than 1 lag, then
%      many of the noisy pixels will be picked out as having significant
%      autocorrelations
%  2) clustering algorithms
%  3) remove non-stationary background 

if nargin < 2
    estNeuronRadius = 3;
    maxNeurons = 100;
elseif nargin < 3
    maxNeurons = 100;
end

estNeuronRadius = round(estNeuronRadius);
estNeuronArea = pi*estNeuronRadius*estNeuronRadius;

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
    fltVideo(:,:,ii) = wiener2(temp,[estNeuronRadius,estNeuronRadius]);
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
        rowVec = max(1, ii-estNeuronRadius):min(width, ii+estNeuronRadius);
        colVec = max(1, jj-estNeuronRadius):min(height, jj+estNeuronRadius);
%         rowVec(rowVec<0) = 0;rowVec(rowVec>width) = 0;rowVec = rowVec(rowVec~=0);
%         colVec(colVec<0) = 0;colVec(colVec<height) = 0;colVec = colVec(colVec~=0);
        for kk=rowVec
            for ll=colVec
                summedCrossCorr(kk,ll) = summedCrossCorr(kk,ll)+...
                    max(xcorr(squeeze(fltVideo(ii,jj,:)),squeeze(fltVideo(kk,ll,:)),maxlag,'coeff'));
                divisor(kk,ll) = divisor(kk,ll)+1;
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

forHistogram = summedCrossCorr(:);
forHistogram = forHistogram(forHistogram ~= 0);
figure();histogram(forHistogram);hold on;
[N,~] = histcounts(forHistogram);
plot(ones(100,1).*threshold,linspace(0,max(N),100),'r','LineWidth',2);
xlabel('Spatial Average Maximum Cross-Correlation');ylabel('Count');
title('Histogram of Maximum Cross-Correlation Coefficients');
legend('Histogram','Bonferroni-Corrected Threshold');

% morphological opening and closing
se = strel('disk',round(estNeuronRadius/4));
se2 = strel('disk',round(estNeuronRadius));
tempBinaryImage = imclose(imopen(binaryCrossCorr,se),se2);

figure();imagesc(tempBinaryImage);colormap('bone');
title('Binary Mask for ROI Detection');

% remove noisy background
maskedVideo = zeros(size(fltVideo));
for ii=1:numFrames 
    maskedVideo(:,:,ii) = fltVideo(:,:,ii).*tempBinaryImage;
end

% h = fspecial('laplacian');
% figure();imagesc(filter2(h,summedCrossCorr,'same'));


cutoff = .2;
%{
% full version of cross-correlation adjacency matrix
% get all possible cross-correlations between nearby pixels
%  linear indexing to row-column indexing
bigMat = zeros(width*height,width*height);
for ii=1:width*height
    [rowInd1,colInd1] = ind2sub(size(tempBinaryImage),ii);
    if tempBinaryImage(rowInd1,colInd1) ~= 0
        rowVec = rowInd1-estNeuronSize:rowInd1+estNeuronSize;
        rowVec(rowVec<0) = 0;rowVec(rowVec>width) = 0;
        rowVec = rowVec(rowVec~=0);
        colVec = colInd1-estNeuronSize:colInd1+estNeuronSize;
        colVec(colVec<0) = 0;colVec(colVec>height) = 0;
        colVec = colVec(colVec~=0);
        for jj=rowVec
            for kk=colVec
                rowInd2 = jj;colInd2 = kk;
                secondInd = sub2ind(size(tempBinaryImage),rowInd2,colInd2);
                bigMat(ii,secondInd) = max(xcorr(squeeze(fltVideo(rowInd1,colInd1,:)),...
                    squeeze(fltVideo(rowInd2,colInd2,:)),maxlag,'coeff'));
            end
        end
    end
end
dissimilarity = 1 - bigMat(tril(true(size(bigMat)), -1)); % lower part of the matrix
dissimilarity = dissimilarity(:);

% clustering works on dissimilarity
Z = linkage(dissimilarity', 'complete');
t = cluster(Z, 'cutoff', cutoff, 'criterion', 'distance');
% t = cluster(Z, 'maxclust', maxNeurons);

col = 'rbmg';
figure(); hold on;
imagesc(tempBinaryImage);
colormap gray;
for i = 1:max(t)
    indexArray = (t == i);
    [r, c] = ind2sub([width, height], find(indexArray));
%     K = boundary(c, r);
%     plot(c(K), r(K), 'linewidth', 2);
   plot(c, r, '.', 'markersize', 15);
end
set(gca, 'YDir', 'reverse');
title('hierarchical clusters (big version)');
%}

%%%% smaller version of adjacency matrix and hierarchical clustering
% store values in a sparse matrix by making an index array and a value
% array
usePixels = find(tempBinaryImage); % only use pixels deemed signal
idx = find(tril(true(length(usePixels)), -1));
[tempr, tempc] = ind2sub([length(usePixels), length(usePixels)], idx);
indexArray = [tempr, tempc];
usePixelArray = usePixels(indexArray);

% distance between pixels
% keep formatting consistent with the output of pdist so that we can use
% linkage and clustering on it
col = ceil(usePixelArray / width);
row = usePixelArray - width * (col - 1);
pixelDist = sqrt((col(:, 1) - col(:, 2)).^2 + (row(:, 1) - row(:, 2)).^2);

% get cross correlations between pairs of pixels
xcorrArray = zeros(length(indexArray), 1);
for ii = 1:length(xcorrArray)
    if pixelDist(ii) > 2*estNeuronRadius
        continue; 
    end
    
    idxr = row(ii, 1);
    idxc = col(ii, 1);
    jdxr = row(ii, 2);
    jdxc = col(ii, 2);
    
    xcorrArray(ii) = max(xcorr(squeeze(fltVideo(idxr, idxc, :)),...
        squeeze(fltVideo(jdxr, jdxc, :)), maxlag, 'coeff'));
end
% clustering works on dissimilarity
dissimilarity = 1 - xcorrArray; % for linkage, smaller means closer together
Z = linkage(dissimilarity', 'complete');
%t = cluster(Z, 'cutoff', cutoff, 'criterion', 'distance');
t = cluster(Z, 'maxclust', maxNeurons);

figure(); hold on;
imagesc(tempBinaryImage);
colormap gray;
for i = 1:max(t)
    indexArray = (t == i);
    [r, c] = ind2sub([width, height], usePixels(indexArray));
%     K = boundary(c, r);
%     plot(c(K), r(K), 'linewidth', 2);
   plot(c, r, '.', 'markersize', 15);
end
set(gca, 'YDir', 'reverse');
title('hierarchical clusters (small version)');

figure(); dendrogram(Z, 0, 'colorthreshold', cutoff);

% bwconncomp finds groups in the binary image
Components = bwconncomp(tempBinaryImage);
Centroids = regionprops(Components,'Centroid');
end
