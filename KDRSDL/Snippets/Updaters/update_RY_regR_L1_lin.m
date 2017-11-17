function [vars] = update_RY_regR_L1_lin(vars, params)
%UPDATE_RY_REGR_L1_LIN Solves for R, L1 penalty
%
% Mehdi Bahri - Imperial College London
% July, 2016
%
% Last modified August, 2017

if params.TIME > 2
    tic
end

S = vars.S;
Uc = vars.A;
Ur = vars.B;
Xt = vars.Xt;
Y = vars.Y;
% Yt = vars.Yt;

mu = params.mu;
% mu_k = vars.mu_k;
r = params.r;
alpha = params.alpha;

% R1 = zeros(r, r, params.Nobs);

% red_S = ttm(tensor(S), {Uc', Ur'}, [1 2]);
% red_S = ( red_S.data + mu_k * vars.K + Yt ) ./ mu_k;

% U = (-mu / mu_k) * (Uc'*Uc);
% V = Ur'*Ur;

% for k=1:params.Nobs
%    R1(:,:,k) = dlyap(U, V, red_S(:,:,k));
% end

% K1 = soft_shrinkage(R1 - Yt / mu_k, alpha/mu_k);

RR = ttm(tensor(vars.R), {Uc'*Uc, Ur'*Ur}, [1, 2]);
Delta = ttm(tensor(vars.Xt + vars.Y ./ mu), {Uc', Ur'}, [1, 2]);

if params.DEGREE_3_REG
    shrink_cst = alpha * norm(vars.A, 'fro') * norm(vars.B, 'fro') /mu;
else
    shrink_cst = alpha / mu;
end

R1 = soft_shrinkage(vars.R - (RR.data - Delta.data) / ( norm(Uc)^2 * norm(Ur)^2 ), shrink_cst);

vars.R = R1;
% vars.K = K1;

L = ttm(tensor(R1), {Uc, Ur}, [1 2]);
L = L.data;

vars.Y = Y + mu*(Xt - L);
% vars.Yt = Yt + mu_k*(K1 - R1);
% vars.mu_k = params.rho * vars.mu_k;

if params.TIME > 2
    up_t_time = toc;
    fprintf('Updated L1 regularized R in %fs (TT)\n', up_t_time);
end

end

