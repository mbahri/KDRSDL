function [vars] = update_B_Fro_split(vars, params)
%COMP_UR_MMX_GROUP_LASSO Optimal Ur with Group Lasso penalty, MMX
%
% Mehdi Bahri - Imperial College London
% June, 2016

% S = vars.S;
% T = vars.T;
% Uc = vars.Uc;
% Y_r = vars.Y_r;
% 
% Urn = mmx('mult', S, Uc, 'tn');
% Urn = sum(mmx('mult', Urn, T) ,3) + vars.mu_r*vars.B + Y_r;
% 
% V = params.mu*(Uc'*Uc);
% 
% Urd = mmx('mult', T, V, 'tn');
% Urd = sum(mmx('mult', Urd, T), 3);
% 
% Ur = Urn / (vars.mu_r * eye(params.r) + Urd);
% vars.Ur = Ur;
% vars.B = solve_l1l2(Ur - Y_r ./ vars.mu_r, params.alpha_r / vars.mu_r);
% 
% vars.Y_r = Y_r + vars.mu_r * (vars.B - Ur); % Check sign
% vars.mu_r = params.rho*vars.mu_r;

UU = vars.U'*vars.U;
Sum_kuuk = zeros(params.r, params.r);
for k=1:params.Nobs
    Sum_kuuk = Sum_kuuk + vars.R(:,:,k)'*UU*vars.R(:,:,k);
end

Red_S = zeros(size(vars.B));
for k=1:params.Nobs
    Red_S = Red_S + vars.S(:,:,k)'*vars.U*vars.R(:,:,k);
end

vars.V = dlyap(...
    eye(size(vars.V, 1)), ...
    -params.mu / vars.mu_v * Sum_kuuk, ...
    vars.B + vars.Y_v / vars.mu_v + Red_S / vars.mu_v ...
);

vars.B = prox_fro(vars.V - vars.Y_v / vars.mu_v, ...
    params.alpha * norm(vars.A, 'fro') * l1_norm(vars.K) / vars.mu_v ...
);

if params.TIME > 2
    fprintf('B updated split %f\n', norm(vars.B, 'fro'));
end

end