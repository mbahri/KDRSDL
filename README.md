# KDRSDL
Code for [Robust Kronecker-Decomposable Component Analysis for Low-Rank Modeling](https://arxiv.org/abs/1703.07886) (ICCV 2017).

Author: [Mehdi Bahri](http://bahri.io)

This code has been adapted from my MSc thesis work at Imperial College (summer of 2016).

## Requirements

* Matlab (tested on Matlab 2015a and later)
* Image Processing Toolbox for noise generation (`imnoise`)
* Control System Toolbox for discrete-time Silverster's equation solver (`dlyap`)

## How to use the code

The main function is kdrsdl.
```matlab
function [A,B, R, L, E] = kdrsdl(foo)
```

## Background subtraction experiments

For background subtraction experiments, the implementation differs slightly from the one presented in the paper in that the mean value of the low-rank tensor's frontal slices is estimated separately and robustly as in [1]. That is, we write ![L = M + L'](http://quicklatex.com/cache3/4b/ql_7ad77138d8fb1f2eea471eb16f7d514b_l3.png) and, at step ![t](http://quicklatex.com/cache3/8d/ql_5121a27906c28c9080bdc88d7480e28d_l3.png), we compute the estimator:

![M_hat^t = (1/N)*Sum_n(X_n - E_n^t)](http://latex.codecogs.com/gif.latex?%5Cmathbf%7B%5Chat%7BM%7D%5Et%7D%20%3D%20%5Cfrac%7B1%7D%7BN%7D%20%5Csum_n%20%28%5Cmathbf%7BX%7D_n%20-%20%5Cmathbf%7BE%7D_n%5Et%29)

We then subtract the current estimated mean, and the outliers, to the input tensor. Line 5 of Algorithm 1 becomes:

![X_tilde^t = X - E^t - M_hat^t](http://quicklatex.com/cache3/31/ql_cbe81c6fee4215c95790a8ea68cef331_l3.png)

All other steps remain identical.

[1] G. Mateos and G. B. Giannakis. Robust PCA as Bilinear 937 Decomposition With Outlier-Sparsity Regularization. IEEE 938 Transactions on Signal Processing, 60(10):5176–5190, 10 939 2012
