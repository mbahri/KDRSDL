function [vars] = init_via_svd(vars, params)
%INIT_VIA_SVD Initialization of the factors by SVD
%
% Mehdi Bahri - Imperial College London
% April, 2016
%
% Last modified August, 2017

if params.TIME > 1
    tic
end

vars.R = zeros(params.r, params.r, params.Nobs);

% Left basis
vars.A = zeros(params.n, params.r);
% Right basis
vars.B = zeros(params.m, params.r);

for i=1:params.Nobs
    [U, S, V] = svd(vars.X(:,:,i));
    vars.A = vars.A + U(:,1:params.r);
    vars.B = vars.B + V(:,1:params.r);
    vars.R(:,:,i) = S(1:params.r, 1:params.r);
end

vars.A = vars.A / params.Nobs;
vars.B = vars.B / params.Nobs;

if params.TIME > 1
    svdtime = toc;
    fprintf('Initialized A, B, R via SVD (%fs)\n', svdtime);
end

end