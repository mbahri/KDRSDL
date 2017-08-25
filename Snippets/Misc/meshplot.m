function [ PP ] = meshplot(M, title_str)
%Mesh plot the input
% saves the plot to a file with name 'fname'

PP = mesh(M);
title(title_str);
view(7, 14);

end

