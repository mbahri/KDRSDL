function [ Data, Info ] = kdrsdl_mv_sub( X, MV, varargin)
%KDRSDL_MV KDRSDL with missing values
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
params.update_E = @(vars, params) (...
    update_E_mv(vars, params)...
);

% params.update_RY = @(vars, params) (...
%     update_RY_regR_L1_lin(vars, params) ...
% );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Additional variables and overrides
vars.K = vars.R;
vars.Yt = zeros(params.r, params.r, params.Nobs);
vars.Omega_bar = boolean(MV);
vars.Omega = boolean(ones(size(MV)) - MV);   % Mask of observed values
params.MV = true;
params.DEGREE_3_REG = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Additional penalty parameters and overrides
vars.mu_k = init_mu_with_norms(vars.R, 1.25, params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Algorithm and specific output values
[Out, Info] = kdrsdl_core(vars, params);

% Assume zero where no observation was made
Out.E(vars.Omega_bar) = 0;

Data.A = Out.A;
Data.B = Out.B;
Data.R = Out.K; % Replacing by Out.R doesn't change much, but the soft-thresholding operator is applied to K
Data.E = Out.E;

if params.MEAN
    Data.M = Out.M(:,:,1);
    Data.L = wm_make_L(Data.A, Data.B, Data.R) + Data.M;
else
    Data.L = wm_make_L(Data.A, Data.B, Data.R);
end

end