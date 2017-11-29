function [vars] = update_RY_regR_L1_d3_sub(vars, params)
%UPDATE_RY_REGR_L1_PRO Solves for R, L1 penalty
%
% Mehdi Bahri - Imperial College London
% November, 2017

if params.TIME > 2
    tic
end

S = vars.S;
Uc = vars.A;
Ur = vars.B;
Xt = vars.Xt;
Y = vars.Y;
Yt = vars.Yt;

mu = params.mu;
mu_k = vars.mu_k;
r = params.r;
alpha = params.alpha;

R1 = zeros(r, r, params.Nobs);

red_S = ttm(tensor(S), {Uc', Ur'}, [1 2]);
red_S = ( red_S.data + mu_k * vars.K + Yt ) ./ mu_k;

U = (-mu / mu_k) * (Uc'*Uc);
V = Ur'*Ur;

for k=1:params.Nobs
   R1(:,:,k) = dlyap(U, V, red_S(:,:,k));
end

K1 = soft_shrinkage(R1 - Yt / mu_k, ... 
    alpha * norm(vars.A, 'fro') * norm(vars.B, 'fro') /mu_k);

vars.R = R1;
vars.K = K1;

L = ttm(tensor(R1), {Uc, Ur}, [1 2]);
L = L.data;

vars.Y = Y + mu*(Xt - L);
vars.Yt = Yt + mu_k*(K1 - R1);
vars.mu_k = params.rho * vars.mu_k;

if params.TIME > 2
    up_t_time = toc;
    fprintf('Updated L1 regularized R in %fs (TT) with deg 3 reg\n', up_t_time);
end

end
