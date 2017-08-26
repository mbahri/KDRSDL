function [ B ] = normalize_dynamic_range( A )
%NORMALIZE_DYNAMIC_RANGE Rescales the image so that the dynamic range is [0
%1]
%
% Mehdi Bahri - Imperial College London
% July, 2016

B = rescale(double(A),0,255);
B = rescale(B,0,1);

end

