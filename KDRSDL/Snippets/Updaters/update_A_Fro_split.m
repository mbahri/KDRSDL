function [vars] = update_A_Fro_split(vars, params)
%UPDATE_A_FRO_SPLIT Compute A with splitting
%
% Mehdi Bahri - Imperial College London
% November, 2017

fprintf('%f\n', norm(vars.A, 'fro'));

% S = vars.S;
% T = vars.T;
% Ur = vars.Ur;
% Y_c = vars.Y_c;

% Ucn = mmx('mult', S, Ur);
% Ucn = sum(mmx('mult', Ucn, T, 'nt') ,3) + vars.mu_c*vars.A + Y_c;
% 
% V = params.mu*(Ur'*Ur);
% 
% Ucd = mmx('mult', T, V);
% Ucd = sum(mmx('mult', Ucd, T, 'nt'), 3);
% 
% Uc = Ucn / (vars.mu_c * eye(params.r) + Ucd);
% vars.Uc = Uc;
% 
% vars.A = solve_l1l2(Uc - Y_c ./ vars.mu_c, params.alpha_c / vars.mu_c);
% 
% vars.Y_c = Y_c + vars.mu_c * (vars.A - Uc); % Need to check the sign
% vars.mu_c = params.rho*vars.mu_c;

VV = vars.V'*vars.V;
Sum_kvvk = zeros(params.r, params.r);
for k=1:params.Nobs
    Sum_kvvk = Sum_kvvk + vars.R(:,:,k)*VV*vars.R(:,:,k)';
end

Red_S = zeros(size(vars.A));
for k=1:params.Nobs
    Red_S = Red_S + vars.S(:,:,k)*vars.V*vars.R(:,:,k)';
end

vars.U = dlyap(...
    eye(size(vars.U, 1)), ...
    -params.mu / vars.mu_u * Sum_kvvk, ...
    vars.A + vars.Y_u / vars.mu_u + Red_S / vars.mu_u ...
);

vars.A = prox_fro(vars.U - vars.Y_u / vars.mu_u, ...
    params.alpha * norm(vars.B, 'fro') * l1_norm(vars.K) / vars.mu_u ...
);

if params.TIME > 2
    fprintf('A updated split %f\n', norm(vars.A, 'fro'));
end

end