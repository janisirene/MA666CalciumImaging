function [SpatMap,CaSignal,Spikes,width,height,Cn,P,options] = CaImSegmentation(VideoFileName,maxNeurons,estNeuronSize)
%CaImSegmentation.m
%   See Pnevmatikakis & Paninski, 2014 & 2016, for their matrix factorization
%     algorithm to automate image segmentation for calcium imaging data.
%     Briefly, if d is the total number of pixels in an image, T is the number
%     of time steps in a video, and K is the maximum number of neurons, then 
%     the observed fluorescence signal Y (matrix size d-by-T) can be
%     regarded as a noisy function of the true signal,F, by Y = F + E, where E is
%     Gaussian noise.  F, in turn, can be regarded as a signal of the form
%     F = AC+B, where A is a d-by-K matrix representing the spatial footprint
%     of neuron k, C is a K-by-T matrix representing the calcium activity of
%     the k-th neuron, and B is a d-by-T matrix representing the background
%     fluorescence at each pixel.  The algorithm infers the matrices
%     A and C from the observed fluorescence signal Y using matrix factorization
%     and convex optimization based on several assumptions:
%       1) the difference between F = AC+B and Y will be the noise E, thus
%        constraining the following steps to find values for A, C, and B
%        that force abs(Y-F) < noise ... they make an estimate of the noise
%        as a first step in the algorithm
%       2) each neuron will be small relative to the size of the original
%        image, of approximate size based on the input estNeuronSize, so 
%        the matrix A should be sparse
%       3) the calcium signal, C, for each neuron will be an AR(p) process
%        of the following form, for an individual neuron, 
%          c(t) = SUM [j=1 to p] gamma-p*c(t-j) + s(t), where s(t) is the
%          number of actual spikes at time t
%
%  This code must be run after downloading ca_source_extraction-master 
%   from the GitHub repository created for the papers referenced above -
%   Go to: https://github.com/epnev/ca_source_extraction.git
%  You also need CVX placed in your MATLAB directory - http://cvxr.com/cvx/download/
%
%INPUT: VideoFileName - file name of the video to be segmented as a string,
%          'Video.tif'
%       maxNeurons - maximum neurons in the video, this is an estimate that
%        constrains the maximum number of neurons, this value will affect
%        the number of neurons found by the algorithm
%       estNeuronSize - estimate of the average neuron diameter in pixels
%
%OUTPUT: SpatMap - d-by-K matrix representing the spatial map for the
%         k-th neuron
%         To view the maps, run the following code:
%          % if d = width*height, in pixels
%          myMap = zeros(width,height,K);
%          for ii=1:K 
%             myMap(:,:,ii) = reshape(SpatMap(:,ii),[width,height]);
%          end
%          % make figure with a representative example
%          figure();imagesc(squeeze(myMap(:,:,1)));colormap jet;
%
%        CaSignal - a K-by-T matrix representing the inferred calcium
%         signal for the k-th neuron across all time steps, dF/F
%        width - width of original image in pixels
%        height - height of original image in pixels
%
% If you are unsure of the neuron size input, you must be very cautious because
% the initialization procedure is heavily dependent on that value.  Check 
% greedyROI2d.m for that procedure, specifically look to imblur()
%Created: 2016/02/09, 24 Cummington, Boston
% Byron Price
%Updated: 2016/11/06
% By: Byron Price

% Run demoMovie with defaults if no arguments given
if nargin == 0 
    VideoFileName = which('demoMovie.tif');
    maxNeurons = 30;
    estNeuronSize = 4;
end
if ischar(VideoFileName)
    [~,~,ext] = fileparts(VideoFileName);
    if strcmpi(ext,'.mat')
        Y = struct2array(load(VideoFileName));
    else
        nam = VideoFileName;          % insert path to tiff stack here
    %     sframe=1;						% user input: first frame to read (optional, default 1)
        Y = readTifStack(nam);
    end
else
    Y = VideoFileName;
end

addpath(genpath('utilities'));
             
Y = Y - min(Y(:)); 
if ~isa(Y,'double');    Y = double(Y);  end         % convert to single

[d1,d2,T] = size(Y);                                % dimensions of dataset
d = d1*d2;                                          % total number of pixels
width = d1;
height = d2;
% Set parameters

K = maxNeurons;                        % number of components to be found - user defined
tau = estNeuronSize;                   % std of gaussian kernel (size of neuron) - 4 is a good start 
p = 2;                                % order of autoregressive system (p = 0 no dynamics, p=1 just decay, p = 2, both rise and decay)
merge_thr = 0.8;                       % merging threshold
tsub = ceil(size(Y,3) / 5000);

options = CNMFSetParms(...                      
    'd1',d1,'d2',d2,...                         % dimensions of datasets
    'search_method','dilate','dist',8,...       % search locations when updating spatial components
    'deconv_method','constrained_foopsi',...    % activity deconvolution method
    'ssub', 1,...                            % spatial downsampling factor (default: 1)
    'tsub', tsub,...                            % temporal downsampling factor (default: 1)    'temporal_iter',2,...                       % number of block-coordinate descent steps 
    'fudge_factor',0.98,...                     % bias correction for AR coefficients
    'merge_thr',merge_thr,...                    % merging threshold
    'maxthr',0.1,...                           % threshold of max value below which values are discarded (default: 0.1)
    'medw',[3,3],...                % size of median filter (default: [3,3])
    'gSig',tau...
    );
% Data pre-processing

[P,Y] = preprocess_data(Y,p);


% fast initialization of spatial components using greedyROI and HALS

[Ain,Cin,bin,fin,center] = initialize_components(Y,K,tau,options,P);  % initialize
refine_components = false;
if refine_components
    [Ain,Cin,center] = manually_refine_components(Y,Ain,Cin,center,Cn,tau,options);
end
% display centers of found components
Cn =  correlation_image(Y); %reshape(P.sn,d1,d2);  %max(Y,[],3); %std(Y,[],3); % image statistic (only for display purposes)
    
% update spatial components
Yr = reshape(Y,d,T);
%clear Y;
[A,b,Cin] = update_spatial_components(Yr,Cin,fin,[Ain,bin],P,options);
if size(A,2) == 0
    error('Zero spatial components found.')
end

% update temporal components
P.p = 0;    % set AR temporarily to zero for speed
[C,f,P,S,YrA] = update_temporal_components(Yr,A,b,Cin,fin,P,options);

% merge found components
[Am,Cm,K_m,merged_ROIs,Pm,Sm] = merge_components(Yr,A,b,C,f,P,S,options);

% repeat
Pm.p = p;    % restore AR value
[A2,b2,Cm] = update_spatial_components(Yr,Cm,f,[Am,b],Pm,options);
if size(A2,2) == 0
    error('Zero spatial components found.')
end
[C2,f2,P2,S2,YrA2] = update_temporal_components(Yr,A2,b2,Cm,f,Pm,options);

% do some plotting
[A_or,C_or,S_or,P_or] = order_ROIs(A2,C2,S2,P2); % order components
K_m = size(C_or,1);
[C_df,~] = extract_DF_F(Yr,A_or,C_or,P_or,options); % extract DF/F values (optional)

SpatMap = A_or(:,1:end-1);
CaSignal = C_df(1:end-1,:);
Spikes = S_or(1:end-1,:);

%contour_threshold = 0.95;                       % amount of energy used for each component to construct contour plot
figure;
[~,~] = plot_contours(A_or,Cn,options,1); % contour plot of spatial footprints
title('Spatiotemporal Correlation Image');

% for noiseless video of calcium dynamics
% myMat = SpatMap*CaSignal;
% maxColor = max(max(CaSignal));
% T = size(CaSignal,2);
% myMat = reshape(myMat,[width,height,T]);
% figure();
% for tt=1:T
%     imagesc(myMat(:,:,tt));colormap jet;colorbar;caxis([0 maxColor]);
%     pause(5/T)
% end
end
