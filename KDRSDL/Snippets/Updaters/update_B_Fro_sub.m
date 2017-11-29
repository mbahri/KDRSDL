function [vars] = update_B_Fro_sub(vars, params)
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

% UU = vars.U'*vars.U;
% Sum_kuuk = zeros(params.r, params.r);
% RR = vars.R;
% for k=1:params.Nobs
%     Sum_kuuk = Sum_kuuk + RR(:,:,k)'*UU*RR(:,:,k);
% end
% 
% Red_S = zeros(size(vars.B));
% SS = vars.S;
% U = vars.U;
% for k=1:params.Nobs
%     Red_S = Red_S + SS(:,:,k)'*U*RR(:,:,k);
% end

UR = reshape( ...
    permute(ttm(tensor(vars.R), {vars.U}, 1), [2, 1, 3]), ...
    [size(vars.R, 2), size(vars.U, 1)*size(vars.R, 3)] ...
);
UR = UR.data;

Sum_kuuk = UR * UR';

Red_S = reshape(permute(vars.S, [2, 1,3]), [size(vars.S, 2), size(vars.S, 1)*size(vars.S, 3)]) * UR';

vars.V = dlyap(...
    eye(size(vars.V, 1)), ...
    -params.mu / vars.mu_v * Sum_kuuk, ...
    vars.B + vars.Y_v / vars.mu_v + Red_S / vars.mu_v ...
);

% vars.V = ( vars.mu_v * vars.B + vars.Y_v + Red_S ) / (vars.mu_v * eye(params.r) + params.mu * Sum_kuuk);

vars.B = prox_fro(vars.V - vars.Y_v / vars.mu_v, ...
    params.alpha * norm(vars.A, 'fro') * l1_norm(vars.K) / vars.mu_v ...
);

vars.Y_v = vars.Y_v + vars.mu_v * (vars.B - vars.V); % Check sign
vars.mu_v = params.rho*vars.mu_v;

if params.TIME > 2
    fprintf('B updated split %f\n', norm(vars.B, 'fro'));
end

end