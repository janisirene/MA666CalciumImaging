function write3DArrayToTif(full, svMovie)
% writes a 3d array [d1 x d2 x time] to a tif
nT = size(full, 3);
for t = 1:nT
    % Note: if you want a bigger simulation that requires saving each time
    % step individually this is the place to do it
    if ~isempty(svMovie)
        imwrite(full(:, :, t) , svMovie,'tif','WriteMode','append');
    end
end
end