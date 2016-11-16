% script to test CaImGetROIs on simulated data
% TO DO:
%   1. systematically test more snr 
%   2. vary threshold in segmentation algorithm to generate full ROC curve
%   3. add metrics to compare ROI footprint shapes
% author: Janis Intoy
% date: November 9, 2016
% modified: Novemeber 10, 2016 - added hit and false alarm rate calcs

%% simulate data
sz = 20; % sz x sz image
dur = 10; % seconds sampled at 30Hz
nROI = 2; % number of cells
snr = 10; % signal to noise ratio
svMovie = ''; % filename to save movie as .tif

[ROI, full] = simulateCalcImg(sz, dur, nROI, snr, svMovie);

%% automatic segmenation algorithm
[finalBinaryImage,Components,Centroids] = CaImGetROIs(full, 5, nROI);

%% plot: compare Components and ROI
truthMap = false(sz, sz);
for i = 1:nROI
    truthMap(ROI(i).indices) = true;
end

figure(); hold on;
imagesc(truthMap);
colormap gray;
for i = 1:Components.NumObjects
    [r, c] = ind2sub([sz, sz], Components.PixelIdxList{i});
    K = boundary(c, r);
    plot(c(K), r(K), 'linewidth', 2);
end
xlim([1, sz]);
ylim([1, sz]);
axis image;
title('detected ROI overlaid on true');
set(gca, 'YDir', 'reverse');

%% was the footprint detected (above threshold?)
X = false(nROI, Components.NumObjects);
for i = 1:nROI
    ind = sub2ind([sz, sz], ROI(i).center(2), ROI(i).center(1));
    for j = 1:Components.NumObjects
        X(i, j) = any(Components.PixelIdxList{j} == ind);
    end
end

% hit and false alarm rate
hit = sum(sum(X, 2) > 0) / nROI;
fa = sum(sum(X, 1) == 0) / Components.NumObjects;
dprime = norminv(hit, 0, 1) - norminv(fa, 0, 1);

% average number of cells contained in detected ROI
avgNPerROI = mean(sum(X, 1)); 

fprintf('hit rate: %1.3f\t\tfalse alarm rate: %1.3f\n', hit, fa);
fprintf('d'' = %1.3f\n', dprime);
fprintf('average number of cells contained in ROI: %1.3f\n', avgNPerROI);