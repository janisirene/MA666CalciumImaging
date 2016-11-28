function n = noise3(sz, pw)
% n = noise3(sz, pw)
% creates a 3d noise matrix of 1/freq^pw (magnitude spectrum)
% sz - dimensions of desired noise (3 dimensions)
% pw - power for noise (1/freq^pw)
% (pw = 0 should make white noise)
% author: Janis Intoy
% November 28, 2016

n = randn(sz); 

im1 = fftshift(fftn(n));
im1 = im1 ./ abs(im1); % magnitude is set to one

[X, Y, Z] = meshgrid(...
    abs((1:sz(2)) - sz(2)/2) + 1,...
    abs((1:sz(1)) - sz(1)/2) + 1,...
    abs((1:sz(3)) - sz(3)/2) + 1);

% create a radial matrix
R = sqrt(X.^2 + Y.^2 + Z.^2);

Fpw = 1 ./ (R.^pw); % 1/frequency^pw (add 1 so that 0 is not included)
im2 = im1 .* Fpw; % magnitudes

n = real(ifftn(ifftshift(im2)));
n = n ./ std(n(:)); % set power equal to 1 for scaling later
end