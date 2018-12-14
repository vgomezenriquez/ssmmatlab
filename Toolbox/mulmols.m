function [beta, M, e] = mulmols(y, Y)
%
% Given the multivariate linear regression model
%
%  y'_t = Y'_t*beta +  epsilon'_t,       t=1,2,...,n,
%
% or, more compactly,
%
%  y = Y*B + E,
%
% this function computes the multiple OLS estimator, its covariance matrix
% and the residuals. The covariance matrix is not multiplied by Sigmar.
%---------------------------------------------------
% USAGE: [beta,M,e]=mulmols(y,Y)
% where:    y      = an (n x m) matrix of y-vectors
%           Y      = matrix of input variables (n x nY)
%---------------------------------------------------
% RETURNS: beta    = an (nY x m) matrix of regression coefficients
%             M    = an (nY x nY) matrix containing (Y'*Y)^{-1}
%             e    = an ((n-nY) x m) matrix containing the residuals
%---------------------------------------------------
%
% Copyright (c) 21 July 2003 by Victor Gomez
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
[n, m] = size(y);
[mY, nY] = size(Y);
if n ~= mY
    error('y and Y must have the same number of rows in mulmols');
end
[Q, R] = qr(Y);
qy = Q' * y;
beta = R(1:nY, :) \ qy(1:nY, :);
M = pinv(R(1:nY, :));
M = M * M';
e = qy(nY+1:n, :);
