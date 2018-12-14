function [Theta, Omega] = ssmspectfac(Gamma)
%
%
%        This function computes the spectral factorization
%        of the covariance generating function
%        G(z) = Gamma(0)+\sum_{i=1}^q Gamma(i)z^{i} + \sum_{i=1}^q Gamma'(i)z^{-i}.
%        That is, it finds a matrix polynomial \Theta(z) = I+\Theta_1z+ ... +\Theta_qz^{q}
%        such that G(z) = \Theta(z)\Omega\Theta'(z^{-1}). This is achieved by solving
%        the Riccati equation for the Kalman filter based on covariance
%        data only
%
%        Sigma = FSigmaF' - (FSigmaH' + G)(HSigmaH' + Gamma(0))^{-1}(FSigmaH' + G)'
%
%        This equation is the resulta of replacing P with -Sigma in the
%        usual DARE equation
%
%        P = FPF' + GQG' - (FPH' + GS)(R + HPH')^{1}(FPH' + GS)',
%
%        given that in this case Pi = P  + Sigma = 0. See Gomez (2016, Sec.
%        5.6)
%
% Input parameters:
% Gamma   = an nxnxm matrix containing the covariances
%           Gamma(0),...,Gamma(q)
% Output parameters:
% K       = the moving average coefficients
% Omega   = the covariance matrix of the innovations
%
% Copyright (c) 21 July 2003 by Victor Gomez
% Ministerio de Hacienda, Direccion Gral. de Presupuestos,
% Subdireccion Gral. de Analisis y P.E.,
% Alberto Alcocer 2, 1-P, D-34, 28046, Madrid, SPAIN.
% Phone : +34-915835439
% E-mail: vgomez@sepg.minhac.es
%
% The author assumes no responsibility for errors or damage resulting from the use
% of this code. Usage of this code in applications and/or alterations of it should
% be referenced. This code may be redistributed if nothing has been added or
% removed and no money is charged. Positive or negative feedback would be appreciated.
%

[~, n, m] = size(Gamma);
q = m - 1;
F = zeros(n*q, n*q);
F(1:(q - 1)*n, n+1:n*q) = eye((q - 1)*n);
H = zeros(n, n*q);
H(1:n, 1:n) = eye(n);
G = zeros(n*q, n);
for i = 1:q
    G((i - 1)*n+1:i*n, :) = Gamma(:, :, i+1);
end
R = Gamma(:, :, 1);
S = eye(n);
Q = zeros(n, n);
J = eye(n);
[~, Theta, Omega] = ss2if(F, G, H, J, Q, S, R);
