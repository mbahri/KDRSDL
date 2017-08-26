function [ O, X ] = yale_sp( scale, level )
%YALE_SP Loads the Yale database, rescales it, and adds S&P noise
%
% Mehdi Bahri - Imperial College London
% July, 2016

load('yaleb10_full_res.mat');
O = X;
clear X cids colidx rowidx;

if ~(scale == 1)
    O = imresize(O, scale);
end

rng('default');
s = rng;
rng(12031992);
X = imnoise(O, 'salt & pepper', level);
rng(s);

end

