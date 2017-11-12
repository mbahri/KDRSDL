function [vars] = comp_Ur_mmx_group_lasso(vars, params)
%COMP_UR_MMX_GROUP_LASSO Optimal Ur with Group Lasso penalty, MMX
%
% Mehdi Bahri - Imperial College London
% June, 2016

S = vars.S;
T = vars.T;
Uc = vars.Uc;
Y_r = vars.Y_r;

Urn = mmx('mult', S, Uc, 'tn');
Urn = sum(mmx('mult', Urn, T) ,3) + vars.mu_r*vars.B + Y_r;

V = params.mu*(Uc'*Uc);

Urd = mmx('mult', T, V, 'tn');
Urd = sum(mmx('mult', Urd, T), 3);

Ur = Urn / (vars.mu_r * eye(params.r) + Urd);
vars.Ur = Ur;
vars.B = solve_l1l2(Ur - Y_r ./ vars.mu_r, params.alpha_r / vars.mu_r);

vars.Y_r = Y_r + vars.mu_r * (vars.B - Ur); % Check sign
vars.mu_r = params.rho*vars.mu_r;

end