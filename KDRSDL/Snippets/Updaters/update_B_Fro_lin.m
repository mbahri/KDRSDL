function [ vars ] = update_B_Fro_lin( vars, params )
%UPDATE_B_FRO_LIN Computes the optimal B with squared Frob. penalty
%
% Mehdi Bahri - Imperial College London
% April, 2016
%
% Last modified August, 2017

S = vars.S;
T = vars.R;
Uc = vars.A;

mu = params.mu;

% Urn = zeros(params.m, params.r);
% Urd = zeros(params.r, params.r);
% 
% if params.PARALLEL
%     parfor k=1:params.Nobs
%         Urn = Urn + S(:,:,k)' * Uc * T(:,:,k);
%         Urd = Urd + mu*T(:,:,k)'*(Uc'*Uc)*T(:,:,k);
%     end
% else
%     for k=1:params.Nobs
%         Urn = Urn + S(:,:,k)' * Uc * T(:,:,k);
%         Urd = Urd + mu*T(:,:,k)'*(Uc'*Uc)*T(:,:,k);
%     end
% end
% 
% vars.B = Urn / (params.alpha_b*eye(params.r) + Urd);

% Compute Lipschitz constant
AA = Uc'*Uc;
Lsum = zeros(params.r, params.r);
for k=1:params.Nobs
    Lsum = Lsum + T(:,:,k)'*AA*T(:,:,k);
end
L_B = 1.01*norm(Lsum, 'fro');

% Compute sum of gradients
Gradients = zeros(size(vars.B));
for k=1:params.Nobs
    Gradients = Gradients + ( (vars.B * T(:,:,k)'*Uc' - S(:,:,k)' / mu) ...
        * Uc * T(:,:,k) );
end

vars.B = prox_fro(vars.B - Gradients / L_B, ...
    ( params.alpha * norm(Uc, 'fro') * l1_norm(vars.K) ) / (mu * L_B) ...
);

if params.TIME > 2
    fprintf('B updated LADMM %f\n', norm(vars.B, 'fro'));
end

end

