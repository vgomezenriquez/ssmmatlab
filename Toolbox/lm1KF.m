function L = lm1KF(x, nd, s, S, p, ps, q, qs, qS)
%
% this function computes the L matrix for outlier detection
% using the fast Kalman filter algorithm. L is the inverse of the
% Cholesky factor of the covariance matrix of the data.

%     INPUTS:
%           x:  an array containing the ARIMA parameter values
%          nd:  an integer, the length of the differenced series
%           s:  seasonality
%           S:  second seasonality
%           p:  AR order
%          ps:  order of the AR of order s
%           q:  order of the regular MA
%          qs:  order of the MA of order s (1 at most)
%          qS:  order of the MA of order S (1 at most)
%  OUTPUTS:
%      L     : a lower triangular matrix
%
% Copyright (c) 21 July 2015 by Victor Gomez
% Ministerio de Hacienda, Direccion Gral. de Presupuestos,
% Subdireccion Gral. de Analisis y P.E.,
% Alberto Alcocer 2, 1-P, D-34, 28046, Madrid, SPAIN.
% Phone : +34-915835439
% E-mail: VGomez@sepg.minhap.es
%
% The author assumes no responsibility for errors or damage resulting from the use
% of this code. Usage of this code in applications and/or alterations of it should
% be referenced. This code may be redistributed if nothing has been added or
% removed and no money is charged. Positive or negative feedback would be appreciated.
%

[phi, alpha, th] = arimapol(x, s, S, p, ps, 0, 0, 0, q, qs, qS);

% set up ARIMA matrices for the model
[Z, T, H, A, Sigma, Xi] = arimam(phi, alpha, th);

%compute lower triangular matrix Theta such that    \Phi*L = \Theta
[Theta, D] = cTheta(nd, T, Sigma, phi, th);

%compute lower triangular matrix Phi
Phi = eye(nd);
plength = size(phi, 2) - 1;
for i = 1:nd - 1
    nphi = min(plength, nd-i);
    mphi = phi(end-1:-1:end-nphi);
    Phi(i+1:i+nphi, i) = mphi';
end
% tic
%the following computation is fast because MATLAB takes care of the form of
%the matrices when using the \ (mldivide) operator.
L = Theta \ Phi;
L = diag(sqrt(D)) \ L;
% toc
L = sparse(L);
