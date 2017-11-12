function [ Data, Info ] = rpca2d_gl_l1( X, varargin)
%RPCA2D_GL_L1 Solves NO2DRPCA with a L1 penalty on T and Group Lasso
% min (1/2)*(alpha_c*||Uc||_F^2 + alpha_c*||Ur||_F^2 + 
%                   alpha_r*sum(||Tn||_F^2) + sum(labmda_n*||En||_1)
% s.t for all n Xn = UcTnUr^T + En
%
% Mehdi Bahri - Imperial College London
% July, 2016

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Defaults for this algorithm
params.alpha = 1e-2;
% params.alpha_c = 1e-3;
% params.alpha_r = 1e-3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Shared behaviour and default parameter initalization
[vars, params] = parameters_common_init( X, params, varargin{:} );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Update functions overrides
params.update_A = @(vars, params) (...
    update_A_Fro_split(vars, params)...
);
params.update_B = @(vars, params) (...
    update_B_Fro_split(vars, params)...
);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Additional variables and overrides
vars.U = vars.A;
vars.V = vars.B;
vars.K = vars.R;

vars.Y_u = zeros(params.n, params.r);
vars.Y_v = zeros(params.m, params.r);
vars.Yt = zeros(params.r, params.r, params.Nobs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Additional penalty parameters and overrides

% params.mu_c = 1.25 / norm(vars.Uc, 'fro');
% params.mu_r = 1.25 / norm(vars.Ur, 'fro');
vars.mu_t = init_mu_with_norms(vars.T, 1.25, params);

vars.mu_u = 1e-3;
vars.mu_v = 1e-3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Algorithm and specific output values
[Out, Info] = rpca2d_core(vars, params);

Data.A = Out.A;
Data.B = Out.B;
Data.R = Out.K;
Data.E = Out.E;

if params.MEAN
    Data.M = Out.M(:,:,1);
    Data.L = wm_make_L(Out.Uc, Out.Ur, Out.T) + Out.M;
else
    Data.L = wm_make_L(Out.Uc, Out.Ur, Out.T);
end 

%     % Test updating Ur with the previous iteration of Uc (i.e before Uc is
%     % updated)
%     Uc = vars.Uc;
%     B = vars.B;
%     [vars.A, vars.Uc, vars.Y_c, params.mu_c] = comp_Uc_mmx_group_lasso(S, vars.T, vars.Ur, vars.A, vars.Y_c, params);
%     [vars.B, vars.Ur, vars.Y_r, params.mu_r] = comp_Ur_mmx_group_lasso(S, vars.T, Uc, B, vars.Y_r, params);
    
%     % Test enforcing the column sparse Uc, Ur
%     [vars.Uc, vars.A, vars.Y_c, params.mu_c] = comp_Uc_mmx_group_lasso(S, vars.T, vars.B, vars.Uc, vars.Y_c, params);
%     [vars.Ur, vars.B, vars.Y_r, params.mu_r] = comp_Ur_mmx_group_lasso(S, vars.T, vars.A, vars.Ur, vars.Y_r, params);

end