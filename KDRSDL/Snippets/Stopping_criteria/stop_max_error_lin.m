function [ emax ] = stop_max_error_lin( vars, params )
%STOP_MAX_ERROR_LIN Computes the maximum reconstruction error of the slices
%
% Mehdi Bahri - Imperial College London
% April, 2016
%
% Last modified August, 2017

X = vars.X;
Xt = vars.Xt;   % Xt is X - E or X - E - M as needed

e_slice = -inf;

Uc = vars.A;
Ur = vars.B;
T = vars.R;

% Reconstruction error defined as the max reconstruction error of the
% individual slices
for k=1:params.Nobs
    e_slice = get_max_error(X(:,:,k), Xt(:,:,k), T(:,:,k), ...
        Uc, Ur, e_slice);
end

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