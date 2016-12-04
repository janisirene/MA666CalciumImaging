function sccRunTests(inputList, tidx, svDir)

params = struct(...
    'size', 30,...
    'duration', 10,...
    'nROI', [],...
    'snr', [],...
    'noisePW', [],...
    'saveMovie', '',...
    'estNeuronRadius', 5);

if ~exist(svDir, 'dir')
    mkdir(svDir); 
end

fid = fopen(inputList, 'r');
for i = 1:tidx
    trialIdx = fscanf(fid, '%i', 1);
    nROI = fscanf(fid, '%i', 1);
    snr = fscanf(fid, '%i', 1);
    noisePW = fscanf(fid, '%i', 1);
end
fclose(fid);

params.nROI = nROI;
params.snr = snr;
params.noisePW = noisePW;

rng(trialIdx); % for using the same simulated data
[out, truthMap, ROI] = testCaImGetROIsOnSimulation(params);

save(fullfile(svDir, sprintf('out_%i.mat', tidx)), 'out', 'truthMap', 'ROI', 'params');
end