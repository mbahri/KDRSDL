function [ vars ] = update_E(vars, params)
%UPDATE_E Solves for E by soft-thresholding
%
% Mehdi Bahri - Imperial College London
% April, 2016
%
% Last modified August, 2017

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

vars.E = soft_shrinkage(E, lambda/mu);

if params.TIME > 2
    fprintf('E updated\n');
end

end

