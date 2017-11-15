function [ Data, Info ] = kdrsdl_pro_split( X, varargin)
%KDRSDL_PRO_SPLIT Kronecker-Decomposable Sparse Dictionary Learning
% min (1/2)*(alpha_a*||A||_F^2 + alpha_b*||B||_F^2 + 
%                   alpha*sum(||Rn||_1) + sum(labmda_n*||En||_1)
% s.t for all n Xn = ARnB^T + En
%
% Mehdi Bahri - Imperial College London - http://bahri.io/
% July, 2016
%
% Last modified August, 2017
% Accepted for publication at ICCV 2017
% Arxiv: https://arxiv.org/abs/1703.07886

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Defaults for this algorithm
params.alpha = 1e-2;

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
params.update_RY = @(vars, params) (...
    update_RY_regR_L1_pro_lin(vars, params) ...
);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Additional variables and overrides
vars.K = vars.R;
vars.Yt = zeros(params.r, params.r, params.Nobs);

vars.U = vars.A;
vars.V = vars.B;
vars.Y_u = zeros(params.n, params.r);
vars.Y_v = zeros(params.m, params.r);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Additional penalty parameters and overrides
vars.mu_k = init_mu_with_norms(vars.R, 1.25, params);

% vars.mu_u = 1e-3;
% vars.mu_v = 1e-3;

vars.mu_u = 1.25 / norm(vars.U, 'fro');
vars.mu_v = 1.25 / norm(vars.V, 'fro');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Algorithm and specific output values
[Out, Info] = kdrsdl_core(vars, params);

Data.A = Out.A;
Data.B = Out.B;
Data.R = Out.R; % Replacing by Out.R doesn't change much, but the soft-thresholding operator returns K
Data.E = Out.E;

if params.MEAN
    Data.M = Out.M(:,:,1);
    Data.L = wm_make_L(Data.A, Data.B, Data.R) + Data.M;
else
    Data.L = wm_make_L(Data.A, Data.B, Data.R);
end

end