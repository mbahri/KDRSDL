function [vars] = comp_Uc_mmx_nuc(vars, params)
%COMP_UC_MMX_NUC Optimal Uc with Nuclear Norm penalty, MMX
%
% Mehdi Bahri - Imperial College London
% June, 2016

S = vars.S;
T = vars.T;
B = vars.Ur;
Uc = vars.A;
Y_c = vars.Y_c;

An = mmx('mult', S, B);
An = sum(mmx('mult', An, T, 'nt') ,3) + vars.mu_c*Uc + Y_c;

V = params.mu*(B'*B);

Ad = mmx('mult', T, V);
Ad = sum(mmx('mult', Ad, T, 'nt'), 3);

A = An / (vars.mu_c * eye(params.r) + Ad);
vars.Uc = A;
vars.A = sv_shrinkage(A - Y_c ./ vars.mu_c, params.alpha_c / vars.mu_c);

vars.Y_c = Y_c + vars.mu_c * (Uc - A);
vars.mu_c = params.rho*vars.mu_c;

end