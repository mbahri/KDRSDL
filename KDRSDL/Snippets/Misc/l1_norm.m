function [ N ] = l1_norm( X )
%L1_NORM The l1 norm of the input array
% 
% Mehdi Bahri
% Imperial College London
% November 2017

N = norm(X(:), 1);

end

