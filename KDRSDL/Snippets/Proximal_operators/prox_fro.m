function [Y] = prox_fro(X, lambda)
% [Y] = prox_fro(X, lambda)
% Proximal operator of the Frobenius norm
%
% Mehdi Bahri
% Imperial College London
% November 2017

if lambda == 0
    
    Y = X;

else
    X = double(X);
    [U, S, V] = svd(X, 'econ');
    S = diag(prox_l2(diag(S), lambda));
    
    Y = U * S * V';
end