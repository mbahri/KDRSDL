function [ err ] = tensor_relative_error( A, B )
%TENSOR_RELATIVE_ERROR Relative error between tensors A and B in L2-norm
%
% Mehdi Bahri - Imperial College London
% July, 2016

A = tensor(A);
B = tensor(B);

err = norm(A - B) / norm(A);


end

