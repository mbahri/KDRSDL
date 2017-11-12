function [ vars ] = update_E_mv(vars, params)
%UPDATE_E_mv Solves for E with missing values by selective-thresholding
%
% Mehdi Bahri - Imperial College London
% October, 2017

E = zeros(params.n, params.m, params.Nobs);

mu = params.mu;
lambda = params.lambda;

if params.PARALLEL
    parfor k=1:params.Nobs
        E(:,:,k) = vars.X(:,:,k) - vars.A*vars.R(:,:,k)*vars.B' + ...
            (1/mu)*vars.Y(:,:,k);
    end
else
    for k=1:params.Nobs
        E(:,:,k) = vars.X(:,:,k) - vars.A*vars.R(:,:,k)*vars.B' + ...
            (1/mu)*vars.Y(:,:,k);
    end
end

vars.E = selective_shrinkage(E, vars.Omega, vars.Omega_bar, lambda/mu);

if params.TIME > 2
    fprintf('E updated\n');
end

end

