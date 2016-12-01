% batchTest script
% as of Dec 1 this one is outdated - moving batch processing to scc so this
% file might be reappropriate to aggregate the data
clear all;

snrs = [3, 5, 10, 30, 50];
nROIs = 1:7;
nTrials = 17;

% snrs = [3, 5];
% nROIs = 1:2;
% nTrials = 2;

params = struct('size', 30,...
    'duration', 10,...
    'nROI', [],...
    'snr', [],...
    'noisePW', 2,...
    'saveMovie', '',...
    'estNeuronRadius', 5);

%
out(nTrials, length(nROIs), length(snrs)) = struct('hitPixels',[],...
    'faPixels', [], 'dprimePixels', [], 'hitSignal', [],...
    'faSignal', [], 'dprimeSignal', [], 'nDetectedPerTrueROI', [],...
    'hitROI', []);
skipped = [0, 0, 0];
ski = 0;
for i = 1:nTrials
    fprintf('%i of %i trials\n', i, nTrials);
    for r = 1:length(nROIs)
        params.nROI = nROIs(r);
        for s = 1:length(snrs)
            params.snr = snrs(s);
            try
                rng(i);
                out(i, r, s) = testCaImGetROIsOnSimulation(params);
            catch err
                ski = ski + 1;
                skipped(ski, :) = [i, r, s];
            end
            close all;
        end
    end
end
save('out_2.mat', 'out', 'skipped', 'snrs', 'nROIs', 'nTrials');
    %}
load('out_2.mat');

%% 
clear combo;
combo(length(nROIs), length(snrs)) = struct();
flds = fields(out);
for s = 1:length(snrs)
    for r = 1:length(nROIs)
        for fi = 1:length(flds)
            temp = [out(:, r, s).(flds{fi})];
            combo(r, s).(flds{fi}) = [nanmean(temp), ...
                nanstd(temp)/sqrt(length(temp))];
        end
    end
end

%% pixels that get through thresholding
sym = '+o*pxsd^v<>h';
figure(1); clf; hold on;
for r = 1:length(nROIs)
    plot(-1, -1, sym(r)); 
end
for r = 1:length(nROIs)
    hitPixels = vertcat(combo(r, :).hitPixels);
    faPixels = vertcat(combo(r, :).faPixels);
    
    %plot(faPixels(:, 1), hitPixels(:, 1), 'o-');
    surface([faPixels(:, 1), faPixels(:, 1)],...
        [hitPixels(:, 1), hitPixels(:, 1)],...
        zeros(length(snrs), 2),...
        [snrs', snrs'],...
        'facecol', 'no',...
        'edgecol', 'interp',...
        'linewidth', 2,...
        'marker', sym(r),...
        'markersize', 10);
end
xlabel('false alarm rate');
ylabel('hit rate');
title('step 1: identifying pixels with signal');
legend(cellfun(@(x) [num2str(x), ' ROI'], num2cell(nROIs),...
    'UniformOutput', false), 'Location', 'southeast');
xlim([0, 1]);
ylim([.95, 1]);
hcb = colorbar;
title(hcb, 'SNR');
set(hcb, 'Ticks', snrs, 'TickLabels', cellfun(@num2str, num2cell(snrs), ...
    'UniformOutput', false));

%% pixels that get through clustering
figure(2); clf; hold on;
for r = 1:length(nROIs)
    plot(-1, -1, sym(r)); 
end
for r = 1:length(nROIs)
    hitSignal = vertcat(combo(r, :).hitSignal);
    faSignal = vertcat(combo(r, :).faSignal);
    
    % i'd prefer a 2d errorbar but :P
    %plot(faSignal(:, 1), hitSignal(:, 1), 'o-');
    surface([faSignal(:, 1), faSignal(:, 1)],...
        [hitSignal(:, 1), hitSignal(:, 1)],...
        zeros(length(snrs), 2),...
        [snrs', snrs'],...
        'facecol', 'no',...
        'edgecol', 'interp',...
        'linewidth', 2,...
        'marker', sym(r),...
        'markersize', 10);
end
xlabel('false alarm rate');
ylabel('hit rate');
title('step 2: pixels that get clustered');
legend(cellfun(@(x) [num2str(x), ' ROI'], num2cell(nROIs),...
    'UniformOutput', false), 'Location', 'southeast');
xlim([0 .25]);
ylim([.2 1]);
hcb = colorbar;
title(hcb, 'SNR');
set(hcb, 'Ticks', snrs, 'TickLabels', cellfun(@num2str, num2cell(snrs), ...
    'UniformOutput', false));

%% number detected ROI with matching true ROI
figure(3); clf; hold on;
for r = 1:length(nROIs)
    temp = vertcat(combo(r, :).nDetectedPerTrueROI);
    plot(snrs, temp(:, 1), ['-', sym(r)],...
        'markersize', 10, 'linewidth', 2);
end
xlabel('SNR');
ylabel('# detected matched to true ROI');
legend(cellfun(@(x) [num2str(x), ' ROI'], num2cell(nROIs),...
    'UniformOutput', false), 'Location', 'northeast');
set(gca, 'XTick', snrs, 'XTickLabels', cellfun(@num2str, num2cell(snrs),...
    'UniformOutput', false), 'XScale', 'log');
title('ratio of detected ROI to matched true ROI');

%% ROI hit rate vs SNR
figure(4); clf; hold on;
for r = 1:length(nROIs)
    temp = vertcat(combo(r, :).hitROI);
    plot(snrs, temp(:, 1), ['-', sym(r)],...
        'markersize', 10, 'linewidth', 2);
end
xlabel('SNR');
ylabel('ROI hit rate');
legend(cellfun(@(x) [num2str(x), ' ROI'], num2cell(nROIs),...
    'UniformOutput', false), 'Location', 'southeast');
set(gca, 'XTick', snrs, 'XTickLabels', cellfun(@num2str, num2cell(snrs),...
    'UniformOutput', false), 'XScale', 'log');
title('ROI hit rate');