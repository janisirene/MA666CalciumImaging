function [] = AVI_TIF_convert(avi_filename)
%% converts AVI files to tiff stack files
%% Usage: AVI_TIF_convert(avi_filename)
% example: AVI_TIF_convert('simulateCalcImg.avi')
% author: Madhura Baxi
% Date: 8th Nov, 2016
obj = VideoReader(avi_filename);
vid = read(obj);
frames = obj.NumberOfFrames;
for k = 1:frames

  imwrite(vid(:,:,:,k),'simulateCalcImg.tif','tif','WriteMode','append');

end
