function [E] = selective_shrinkage(X, Omega, Omega_bar, lambda)
%SELECTIVE_SHRINKAGE Elementwise soft-shrinkage operator
%
% Apply soft thresholding to the observed values, otherwise keep the
% reconstructed values.
%
% Mehdi Bahri - Imperial College London
% October, 2017

% ones(size(X)) - Omega gives the boolean complement of Omega

X_obs = X .* Omega;                             % Omega(X) observed
X_unk = X .* Omega_bar;                         % Omega_bar(X)

E = ( max(abs(X_obs) - lambda, 0) .* sign(X_obs) ) + X_unk;

end