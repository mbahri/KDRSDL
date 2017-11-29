function [vars] = update_A_Fro_sub(vars, params)
%UPDATE_A_FRO_SPLIT Compute A with splitting
%
% Mehdi Bahri - Imperial College London
% November, 2017

% fprintf('%f\n', norm(vars.A, 'fro'));

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

% VV = vars.V'*vars.V;
% Sum_kvvk = zeros(params.r, params.r);
% RR = vars.R;
% for k=1:params.Nobs
%     Sum_kvvk = Sum_kvvk + RR(:,:,k)*VV*RR(:,:,k)';
% end
% 
% Red_S = zeros(size(vars.A));
% SS = vars.S;
% V = vars.V;
% for k=1:params.Nobs
%     Red_S = Red_S + SS(:,:,k)*V*RR(:,:,k)';
% end

RV = reshape( ...
        ttm(tensor(vars.R), {vars.V}, 2), ...
        [size(vars.R, 2), size(vars.V, 1) * size(vars.R, 3)] ...
);
RV = RV.data;

Sum_kvvk = RV * RV';

Red_S = reshape(vars.S, [size(vars.S, 1), size(vars.S, 2)*size(vars.S, 3)]) * RV';

vars.U = dlyap(...
    eye(size(vars.U, 1)), ...
    -params.mu / vars.mu_u * Sum_kvvk, ...
    vars.A + vars.Y_u / vars.mu_u + Red_S / vars.mu_u ...
);

% vars.U = (vars.mu_u * vars.A + vars.Y_u + Red_S) / (vars.mu_u * eye(params.r) + params.mu * Sum_kvvk);

vars.A = prox_fro(vars.U - vars.Y_u / vars.mu_u, ...
    params.alpha * norm(vars.B, 'fro') * l1_norm(vars.K) / vars.mu_u ...
);

vars.Y_u = vars.Y_u + vars.mu_u * (vars.A - vars.U); % Need to check the sign
vars.mu_u = params.rho*vars.mu_u;

if params.TIME > 2
    fprintf('A updated split %f\n', norm(vars.A, 'fro'));
end

end