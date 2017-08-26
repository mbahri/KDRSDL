%LOAD_HALL Loads the hall test video and rescales the dynamic range
%
% Mehdi Bahri - Imperial College London
% July, 2016

load('hall.mat');
O = normalize_dynamic_range(vid2);
GT = normalize_dynamic_range(GT);
GT_frames = GT_frames - frames(1) + 1;
clear vid1 vid2 frames imHeight imWidth nFrames