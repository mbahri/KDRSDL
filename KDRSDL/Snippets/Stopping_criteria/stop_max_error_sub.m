function [ emax ] = stop_max_error_sub( vars, params )
%STOP_MAX_ERROR_SUB Computes the maximum reconstruction error of the slices
%
% Mehdi Bahri - Imperial College London
% April, 2016
%
% Last modified August, 2017

X = vars.X;
Xt = vars.Xt;   % Xt is X - E or X - E - M as needed

e_slice = -inf;
e_core = -inf;
e_bases = -inf;

Uc = vars.A;
Ur = vars.B;
T = vars.K;

% If degree 3 and sub then there are auxiliary variables for A and B
e_A = matrix_relative_error(vars.A, vars.U);
e_B = matrix_relative_error(vars.B, vars.V);
if params.DEGREE_3_REG
    e_bases = max(matrix_relative_error(vars.A, vars.U), matrix_relative_error(vars.B, vars.V));
end

% Reconstruction error defined as the max reconstruction error of the
% individual slices
for k=1:params.Nobs
    e_slice = get_max_error(X(:,:,k), Xt(:,:,k), T(:,:,k), ...
        Uc, Ur, e_slice);
end

% Splitting error on the core
for k=1:params.Nobs
    e_core = max(e_core, matrix_relative_error(vars.K(:,:,k), vars.R(:,:,k)));
end

fprintf('e_slice = %f | e_core = %f | e_A = %f | e_B = %f\n', e_slice, e_core, e_A, e_B);
emax = max([e_slice, e_core, e_bases]);

end

function [emax] = get_max_error(X, Xt, T, Uc, Ur, currmax)

Error = norm(Xt - Uc*T*Ur', 'fro') /...
    norm(X, 'fro');

if Error > currmax
   emax = Error;
else
   emax = currmax;
end

end