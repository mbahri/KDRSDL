function [ vars ] = update_A_Fro_lin( vars, params )
%UPDATE_A_FRO_LIN Computes the optimal A with Frob. penalty by LADMM
%
% Mehdi Bahri - Imperial College London
% April, 2016
%
% Last modified August, 2017

% fprintf('%f\n', norm(vars.A, 'fro'));

S = vars.S;
T = vars.R;
Ur = vars.B;

mu = params.mu;

% Ucn = zeros(params.n, params.r);
% Ucd = zeros(params.r, params.r);
% 
% if params.PARALLEL
%     parfor k=1:params.Nobs
%         Ucn = Ucn + S(:,:,k) * Ur * T(:,:,k)';
%         Ucd = Ucd + mu*T(:,:,k)*(Ur'*Ur)*T(:,:,k)';
%     end
% else
%     for k=1:params.Nobs
%         Ucn = Ucn + S(:,:,k) * Ur * T(:,:,k)';
%         Ucd = Ucd + mu*T(:,:,k)*(Ur'*Ur)*T(:,:,k)';
%     end
% end
%     
% vars.A = Ucn / (params.alpha_a*eye(params.r) + Ucd);

% Compute Lipschitz constant
BB = Ur'*Ur;
Lsum = zeros(params.r, params.r);
for k=1:params.Nobs
    Lsum = Lsum + T(:,:,k)*BB*T(:,:,k)';
end
L_A = 1.01*norm(Lsum, 'fro');

% Compute sum of gradients
Gradients = zeros(size(vars.A));
for k=1:params.Nobs
    Gradients = Gradients + ( (vars.A * T(:,:,k)*Ur' - S(:,:,k) / mu) ...
        * Ur * T(:,:,k)' );
end

vars.A = prox_fro(vars.A - Gradients / L_A, ...
    ( params.alpha * norm(Ur, 'fro') * l1_norm(vars.K) ) / (mu * L_A) ...
);

if params.TIME > 2
    fprintf('A updated LADMM %f\n', norm(vars.A, 'fro'));
end
    
end

