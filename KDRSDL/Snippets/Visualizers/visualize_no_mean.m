function [ ] = visualize_no_mean( vars, params )
%VISUALIZE_NO_MEAN Visualize the progress of the training process
%
% Mehdi Bahri - Imperial College London
% May, 2016
%
% Last modified August, 2017

Uc = vars.A;
Ur = vars.B;
T = vars.R;

subplot(2,3,1), imshow(vars.X(:,:,1), []), title('First image')
subplot(2,3,2), imshow(Uc*T(:,:,1)*Ur', []), title('First reconstruction')
subplot(2,3,3), imshow(vars.E(:,:,1), []), title('First outliers')

subplot(2,3,4), stemplot(Uc, 'Spectrum of A');
subplot(2,3,5), stemplot(Ur, 'Spectrum of B');
subplot(2,3,6), meshplot(T(:,:,1), 'Sparsity of first representation');
drawnow;
    
end

