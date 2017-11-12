function [ y ] = prox_l2( x, lambda )
%PROX_L2 Proximal operator of the l2 (not squared) norm
%
% Mehdi Bahri
% Imperial College London
% November 2017

n = norm(x);
y = max(n - lambda, 0) * x / n;

end

