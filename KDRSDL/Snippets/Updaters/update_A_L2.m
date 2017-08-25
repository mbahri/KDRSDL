function [ vars ] = update_A_L2( vars, params )
%UPDATE_A_L2 Computes the optimal A with squared Frob. penalty
%
% Mehdi Bahri - Imperial College London
% April, 2016
%
% Last modified August, 2017

S = vars.S;
T = vars.R;
Ur = vars.B;

mu = params.mu;

Ucn = zeros(params.n, params.r);
Ucd = zeros(params.r, params.r);

if params.PARALLEL
    parfor k=1:params.Nobs
        Ucn = Ucn + S(:,:,k) * Ur * T(:,:,k)';
        Ucd = Ucd + mu*T(:,:,k)*(Ur'*Ur)*T(:,:,k)';
    end
else
    for k=1:params.Nobs
        Ucn = Ucn + S(:,:,k) * Ur * T(:,:,k)';
        Ucd = Ucd + mu*T(:,:,k)*(Ur'*Ur)*T(:,:,k)';
    end
end
    
vars.A = Ucn / (params.alpha_a*eye(params.r) + Ucd);

if params.TIME > 2
    fprintf('A updated\n');
end
    
end

