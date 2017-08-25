function [ PP ] = stemplot(M, title_str)
%STEMPLOT Stem plot the singular values of the input
% saves the plot to a file with name 'fname'

[~, S, ~] = svd(M);
S = diag(S);

PP = stem(S);
title(title_str);

end

