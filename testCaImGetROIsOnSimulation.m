% script to test CaImGetROIs on simulated data
% TO DO:
%   1. systematically test a variety of snr 
%   2. add metrics to compare truth and output (ROC curve requires changing
%   threshold in autosegmentation algorithm, perhaps just hit and FA rates
%   for now?)

%% simulate data
sz = 100; % sz x sz image
dur = 10; % seconds sampled at 30Hz
nROI = 20; % number of cells
snr = 3; % signal to noise ratio
svMovie = ''; % filename to save movie as .avi

[ROI, full] = simulateCalcImg(sz, dur, nROI, snr, svMovie);

%% automatic segmenation algorithm
[finalBinaryImage,Components,Centroids] = CaImGetROIs(full);

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

%% 