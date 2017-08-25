function [ err ] = matrix_relative_error( A, B )
%MATRIX_RELATIVE_ERROR Relative error between matrices A and B in Frob. norm
%
% Mehdi Bahri - Imperial College London
% July, 2016

err = norm(A - B, 'fro') / norm(A, 'fro');

end

