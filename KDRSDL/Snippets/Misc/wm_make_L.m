function [ L ] = wm_make_L( Uc, Ur, T )
%WM_MAKE_L Reconstructs the low-rank component from the bases and the core
%tensor
%
% Mehdi Bahri - Imperial College London
% July, 2016

L = double(ttm(tensor(T), {Uc Ur}, [1 2]));

end

