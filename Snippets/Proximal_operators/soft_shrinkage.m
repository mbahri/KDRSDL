function [E] = soft_shrinkage(X, lambda)
%SOFT_SHRINKAGE Elementwise soft-shrinkage operator
%
% Mehdi Bahri - Imperial College London
% April, 2016

E = max(abs(X) - lambda, 0) .* sign(X);

end