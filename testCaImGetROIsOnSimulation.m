function out = testCaImGetROIsOnSimulation(params)
% TO DO:
%   1. systematically test more snr 
%   2. vary threshold in segmentation algorithm to generate full ROC curve
%   3. add metrics to compare ROI footprint shapes
% author: Janis Intoy
% date: November 9, 2016
% modified: Novemeber 10, 2016 - added hit and false alarm rate calcs

%% parameters
if ~exist('params', 'var')
    sz = 30; % sz x sz image
    dur = 10; % seconds sampled at 30Hz
    nROI = 5; % number of cells
    snr = 30; % signal to noise ratio
    svMovie = ''; % filename to save movie as .tif
    
    estNeuronRadius = 5;
    

else
    sz = params.size;
    dur = params.duration;
    nROI = params.nROI;
    snr = params.snr;
    svMovie = params.saveMovie;
    estNeuronRadius = params.estNeuronRadius;
end

%% simulate data
[ROI, full] = simulateCalcImg(sz, dur, nROI, snr, svMovie);

%% automatic segmenation algorithm
[finalBinaryImage,Components,Centroids, clusterData] = CaImGetROIs(full, estNeuronRadius, nROI);

%% truth map
truthMap = zeros(sz);
for i = 1:nROI
    truthMap(ROI(i).indices) = i; % if overlap the later index dominates
end

%% plot: compare Components and ROI
detectedROI = struct('indices', [], 'trueROI', []);
detectedMap = zeros(sz);
cnt = 1;

figure(); hold on;
imagesc(truthMap);
colormap jet;
for i = 1:max(clusterData.class)
    idx = clusterData.class == i;
    if sum(idx) > 2*estNeuronRadius
        sIdx = clusterData.signalPixels(idx);
        [r, c] = ind2sub([sz, sz], sIdx);
        K = boundary(c, r);
        plot(c(K), r(K), 'k', 'linewidth', 2);
        
        % match with a true ROI
        truth = truthMap(sIdx);

        detectedROI(cnt).trueROI = mode(truth); % will be zero if no match
        detectedROI(cnt).indices = sIdx;
        detectedMap(sIdx) = 1;
        cnt = cnt + 1;
    end
end
% plotTrueROIOverlay(ROI);
axis image;
xlim([1, sz]); ylim([1, sz]);
title('detected ROI overlaid on true');
set(gca, 'YDir', 'reverse');


%% individual pixel analysis (test of noise removal process)
hitPixels = sum(truthMap(:) > 0 & finalBinaryImage(:)) / sum(truthMap(:) > 0);
faPixels = sum(truthMap(:) == 0 & finalBinaryImage(:)) / sum(truthMap(:) == 0);
dprimePixels = norminv(hitPixels, 0, 1) - norminv(faPixels, 0, 1);

%% signal pixels (test of what made it through clustering)
hitSignal = sum(truthMap(:) > 0 & detectedMap(:)) / sum(truthMap(:) > 0);
faSignal = sum(truthMap(:) == 0 & detectedMap(:)) / sum(truthMap(:) == 0);
dprimeSignal = norminv(hitSignal, 0, 1) - norminv(faSignal, 0, 1);

%% ROI analysis
foundROIs = [detectedROI.trueROI];
actualNROI = length(unique(truthMap(:))) - any(truthMap(:) == 0);
nDetectedPerTrue = sum(foundROIs > 0) / actualNROI;

nDetected = length(unique(foundROIs)) - any(foundROIs == 0);
hitROI = nDetected / actualNROI;
% is there a right way to calculate a false alarm rate on ROI? I guess
% that's what the individual pixel way is for 

%% save output variables
out = struct('hitPixels', hitPixels,...
    'faPixels', faPixels,...
    'dprimePixels', dprimePixels,...
    'hitSignal', hitSignal,...
    'faSignal', faSignal,...
    'dprimeSignal', dprimeSignal,...
    'nDetectedPerTrueROI', nDetectedPerTrue,...
    'hitROI', hitROI);
end