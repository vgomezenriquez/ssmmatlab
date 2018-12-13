function [beta, M, e] = bmols(y, Y)
%**************************************************************************
%
% This function computes the OLS estimator, its covariance matrix and the
% white noise residuals. The covariance matrix is not multiplied by
% sigma^2.
%
%   INPUTS:
%       y : data vector
%       Y : matrix with regression variables
%
%  OUTPUTS:
%    beta : OLS estimator
%       M : covariance matrix of beta
%       e : residuals
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
%**************************************************************************

n = length(y);
[mY, nY] = size(Y);
if n ~= mY
    error('y and Y must have the same number of rows in bmols');
end
[Q, R] = qr(Y);
qy = Q' * y;
% beta=R(1:nY,:)\qy(1:nY);
M = pinv(R(1:nY, :));
beta = M * qy(1:nY);
M = M * M';
e = qy(nY+1:n);
