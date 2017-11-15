function [ emax ] = stop_max_error( vars, params )
%STOP_MAX_ERROR Computes the maximum reconstruction error of the slices
%
% Mehdi Bahri - Imperial College London
% April, 2016
%
% Last modified August, 2017

X = vars.X;
Xt = vars.Xt;   % Xt is X - E or X - E - M as needed

e_slice = -inf;
% e_core = -inf;

Uc = vars.A;
Ur = vars.B;
T = vars.R;

% Reconstruction error defined as the max reconstruction error of the
% individual slices
for k=1:params.Nobs
    e_slice = get_max_error(X(:,:,k), Xt(:,:,k), T(:,:,k), ...
        Uc, Ur, e_slice);
end

% Splitting error on the core
% for k=1:params.Nobs
%     e_core = max(e_core, matrix_relative_error(vars.K(:,:,k), vars.R(:,:,k)));
% end

emax = max(e_slice);

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