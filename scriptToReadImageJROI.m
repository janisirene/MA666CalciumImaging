%% read in ImageJ ROI files
% author: Janis Intoy

addpath('ReadImageJROI');

%%%%%% change paths to zip file and tif
ptToROIFiles = fullfile('/Users', 'janisintoy', 'Desktop', ...
    'RoiSet_day15-area3.zip');
ptToTif = fullfile('/Users', 'janisintoy', 'Desktop',...
    'XL011-10_day15_area3_002_1.tif');

sROI = ReadImageJROI(ptToROIFiles);
%vnRectBounds: ['nTop', 'nLeft', 'nBottom', 'nRight']

%% what is the corresponding file for these data?
Y = readTifStack(ptToTif);
Y = double(Y);

%% play tif with ROI?
cmax = max(Y(:));
cmin = min(Y(:));

figure(1); clf; hold on;
h1 = imagesc(Y(:, :, 1));
caxis([cmin, cmax]);
colormap gray;
axis image;

% plot ROIs
for i = 1:length(sROI)
    pos = sROI{i}.vnRectBounds;
    x = [pos(2), pos(2), pos(4), pos(4), pos(2)];
    y = [pos(3), pos(1), pos(1), pos(3), pos(3)];
    
    plot(x, y, 'r', 'linewidth', 1);
end
xlim([1, size(Y, 2)]);
ylim([1, size(Y, 1)]);

for t = 1:size(Y, 3)
    set(h1, 'CData', Y(:, :, t));
    caxis([cmin, cmax]);
    
    pause(.01);
end