function [ vars ] = update_B_L2( vars, params )
%UPDATE_B_L2 Computes the optimal B with squared Frob. penalty
%
% Mehdi Bahri - Imperial College London
% April, 2016
%
% Last modified August, 2017

S = vars.S;
T = vars.R;
Uc = vars.A;

mu = params.mu;

Urn = zeros(params.m, params.r);
Urd = zeros(params.r, params.r);

if params.PARALLEL
    parfor k=1:params.Nobs
        Urn = Urn + S(:,:,k)' * Uc * T(:,:,k);
        Urd = Urd + mu*T(:,:,k)'*(Uc'*Uc)*T(:,:,k);
    end
else
    for k=1:params.Nobs
        Urn = Urn + S(:,:,k)' * Uc * T(:,:,k);
        Urd = Urd + mu*T(:,:,k)'*(Uc'*Uc)*T(:,:,k);
    end
end

vars.B = Urn / (params.alpha_b*eye(params.r) + Urd);

if params.TIME > 2
    fprintf('B updated\n');
end

end

