function [ mu ] = init_mu_with_norms( X, rescaling_coeff, params )
%INIT_MU_WITH_NORMS Dual step initialization
%
% Mehdi Bahri - Imperial College London
% April, 2016

if params.TIME > 1
    tic
end

snorms = 0;
for k=1:params.Nobs
    snorms = snorms + norm(X(:,:,k), 'fro');
end

% Penalty parameter
mu = rescaling_coeff * params.Nobs / snorms; % this one can be tuned

if params.TIME > 1
    normstime = toc;
    fprintf('Initialized mu to %f with norms (%fs)\n', ...
        mu, normstime);
end

end

