# KDRSDL
Code for [Robust Kronecker-Decomposable Component Analysis for Low-Rank Modeling](https://arxiv.org/abs/1703.07886) (ICCV 2017).

Author: [Mehdi Bahri](http://bahri.io)

This code has been adapted from my MSc thesis work at Imperial College (summer of 2016).

## Requirements

* Matlab (tested on Matlab 2015a and later)
* Image Processing Toolbox for noise generation (`imnoise`)
* Control System Toolbox for discrete-time Silverster's equation solver (`dlyap`)
* Optional: Parallel Computing Toolbox

## How to use the code

The main function is `kdrsdl`.
```matlab
function [ Data, Info ] = kdrsdl( X, varargin)
```

`Data` and `Info` are structures that contain, respectively, the recovered components of the decomposition, and some information to monitor convergence (number of iterations and value of the stopping criterion at each step).

The `Experiments` folder contains sample data sets along with helper functions.

### Parameters
Reference of the supported options.
#### 'r' (integer)
Upper-bound on the mode-1 and mode-2 ranks, the maximum dimension of the code matrices `Rn` is `r*r`. Defaults to `min(m, n)` for `m*n` matrices.
#### 'lambda' (double)
The `lambda` parameter in the optimization problem. Defaults to `1 / sqrt(N * max(m, n))` for `m*n*N` tensors.
#### 'tol' (double)
Tolerance for convergence, the algorithm will stop when the stopping criterion is `< tol`. Defaults to `1e-7`.
#### 'maxiter' (integer)
Maximum number of iterations. Defaults to `150`.
#### 'rho' (double)
The multiplicative factor by which the dual-steps `mu` and `mu_k` are increased. Defaults to `1.2`
#### 'enforce_monotonicity' (boolean)
If true, the algorithm will stop as soon as the reconstruction error increases (regardless of weather it would have later decreased). Defaults to `false`
#### 'alpha_a' (double)
Allows giving a weight different from 1 to `||A||_F^2` in the optimization problem. Defaults to `1`.
#### 'alpha_b' (double)
Allows giving a weight different from 1 to `||B||_F^2` in the optimization problem. Defaults to `1`.
#### 'alpha' (double)
The `alpha` parameter in the optimization problem. Defaults to `1e-2`.
#### 'mu' (double)
Allows overriding the initialization of `mu`.
#### 'ground_o' (array of doubles)
Allows feeding in the ground-truth low-rank tensor (useful to write visualizations to monitor training).
#### 'ground_e' (array of doubles)
Allows feeding in the ground-truth sparse tensor (useful to write visualizations to monitor training).
#### 'mean' (boolean)
Robustly estimate the sample mean (see last section of this readme). Defaults to `false`
#### 'parallel' (boolean)
Run some of the updates in parallel with the _Matlab Parallel Computing Toolbox_. Used for computations that are identical and independent on all frontal slices of a tensor. Default to `false` (**warning: enabling it may actually lead to decreased performance due to overhead**).
#### 'time' (boolean)
Level of details with which to display timing information.
* 0: No timing
* 1: Display time taken by each iteration
* 2: Add time spent in intialization to 1
* 3: Add time spent updating the codes `R` to 2, notifies when `A`, `B`, `E` have been updated

### Denoising experiments

We provide both the _Yale-B_ and _Facade_ benchmarks, the procedure is the same in both cases:
```matlab
% Load the original and noisy image(s)
% [Original, Noisy] = yale_sp(rescaling_factor, noise_level)
% Call facade_sp for the Facade benchmark

[O, X] = yale_sp(1, 0.3);

% See parameters section for the supported options
[Data, Info] = kdrsdl(X);

% Display reconstruction
imshow(Data.L(:,:,1), []);
```

### Background subtraction experiments

We provide the _Hall_ data set and two ground-truth frames.
```matlab
% Load the frames and the ground-truth
load_hall
% There are now three variables:
% * O: The frames
% * GT: The ground-truth frames
% * GT_frames: The ids of the two frames that match the ground truth

% See below for explanations about the mean parameter
[Data, Info] = kdrsdl(O, 'mean', 1)
```

## Background subtraction experiments

For background subtraction experiments, the implementation differs slightly from the one presented in the paper in that the mean value of the low-rank tensor's frontal slices is estimated separately and robustly as in [1]. That is, we write ![L = M + L'](http://quicklatex.com/cache3/4b/ql_7ad77138d8fb1f2eea471eb16f7d514b_l3.png) and, at step ![t](http://quicklatex.com/cache3/8d/ql_5121a27906c28c9080bdc88d7480e28d_l3.png), we compute the estimator:

![M_hat^t = (1/N)*Sum_n(X_n - E_n^t)](http://latex.codecogs.com/gif.latex?%5Cmathbf%7B%5Chat%7BM%7D%5Et%7D%20%3D%20%5Cfrac%7B1%7D%7BN%7D%20%5Csum_n%20%28%5Cmathbf%7BX%7D_n%20-%20%5Cmathbf%7BE%7D_n%5Et%29)

We then subtract the current estimated mean, and the outliers, to the input tensor. Line 5 of Algorithm 1 becomes:

![X_tilde^t = X - E^t - M_hat^t](http://quicklatex.com/cache3/31/ql_cbe81c6fee4215c95790a8ea68cef331_l3.png)

All other steps remain identical.

[1] G. Mateos and G. B. Giannakis. Robust PCA as Bilinear 937 Decomposition With Outlier-Sparsity Regularization. IEEE 938 Transactions on Signal Processing, 60(10):5176–5190, 10 939 2012
