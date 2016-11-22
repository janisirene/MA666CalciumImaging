% batchTest script
snrs = [3, 5, 10, 30, 50];
nROIs = 1:7;
nTrials = 50;

snrs = [3, 5];
nROIs = 1:2;
nTrials = 2;

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