function [ r ] = estim_rank( A, t )
%ESTIM_RANK Estimates the rank of A based on the SVD
%
% Estimates the rank by counting the singular values above a given
% threshold t.
%
% Mehdi Bahri - Imperial College London
% July, 2016

[~, S, ~] = svd(A);
S = diag(S);

if nargin < 2
    r = sum(S ./ max(S) > 1e-5);
else
    r = sum(S ./ max(S) > t);
end

end

