function [vars] = comp_Ur_mmx_nuc(vars, params)
%COMP_UR_MMX_NUC Optimal Ur with Group Lasso penalty, MMX
%
% Mehdi Bahri - Imperial College London
% June, 2016

S = vars.S;
T = vars.T;
A = vars.Uc;
Ur = vars.B;
Y_r = vars.Y_r;

Bn = mmx('mult', S, A, 'tn');
Bn = sum(mmx('mult', Bn, T) ,3) + vars.mu_r*Ur + Y_r;

V = params.mu*(A'*A);

Bd = mmx('mult', T, V, 'tn');
Bd = sum(mmx('mult', Bd, T), 3);

B = Bn / (vars.mu_r * eye(params.r) + Bd);
vars.Ur = B;
vars.B = sv_shrinkage(B - Y_r ./ vars.mu_r, params.alpha_r / vars.mu_r);

vars.Y_r = Y_r + vars.mu_r * (Ur - B);
vars.mu_r = params.rho*vars.mu_r;

end