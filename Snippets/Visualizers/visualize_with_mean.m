function [ ] = visualize_with_mean( vars, params )
%VIZUALIZE_WITH_MEAN  Visualize the progress of the training process with
% mean estimation
%
% Mehdi Bahri - Imperial College London
% May, 2016
%
% Last modified August, 2017

Uc = vars.A;
Ur = vars.B;
T = vars.K;

subplot(2,4,1), imshow(vars.X(:,:,1), []), title('First image')
subplot(2,4,2), imshow(Uc*T(:,:,1)*Ur', []), title('First reconstruction minus mean')
subplot(2,4,3), imshow(vars.M(:,:,1), []), title('Estimated mean')
subplot(2,4,4), imshow(vars.E(:,:,1), []), title('First outliers')

subplot(2,4,5); imshow(Uc*T(:,:,1)*Ur' + vars.M(:,:,1), []); title('First reconstruction');
subplot(2,4,6), stemplot(Uc, 'Spectrum of A');
subplot(2,4,7), stemplot(Ur, 'Spectrum of B');
subplot(2,4,8), meshplot(T(:,:,1), 'Sparsity of first representation');
% if isfield(params, 'ground_O')
%     subplot(2,4,8); imshow(params.g_O(:,:,1), []); title('Original')
% end
drawnow;
    
end

