% batchTest script
clear all;

snrs = [3, 5, 10, 30, 50];
nROIs = 1:7;
nTrials = 50;

% snrs = [3, 5];
% nROIs = 1:2;
% nTrials = 2;

params = struct('size', 30,...
    'duration', 10,...
    'nROI', [],...
    'snr', [],...
    'saveMovie', '',...
    'estNeuronRadius', 5);

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
save('out.mat', 'out', 'skipped', 'snrs', 'nROIs', 'nTrials');

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
figure(1); clf; hold on;
for r = 1:length(nROIs)
    hitPixels = vertcat(combo(r, :).hitPixels);
    faPixels = vertcat(combo(r, :).faPixels);
    
    % i'd prefer a 2d errorbar but :P
    plot(faPixels(:, 1), hitPixels(:, 1), 'o-');
end
xlabel('false alarm rate');
ylabel('hit rate');
title('step 1: identifying pixels with signal');
legend(cellfun(@(x) [num2str(x), ' ROI'], num2cell(nROIs),...
    'UniformOutput', false), 'Location', 'southeast');

%% pixels that get through clustering
figure(2); clf; hold on;
for r = 1:length(nROIs)
    hitSignal = vertcat(combo(r, :).hitSignal);
    faSignal = vertcat(combo(r, :).faSignal);
    
    % i'd prefer a 2d errorbar but :P
    plot(faSignal(:, 1), hitSignal(:, 1), 'o-');
end
xlabel('false alarm rate');
ylabel('hit rate');
title('step 2: pixels that get clustered');
legend(cellfun(@(x) [num2str(x), ' ROI'], num2cell(nROIs),...
    'UniformOutput', false), 'Location', 'southeast');