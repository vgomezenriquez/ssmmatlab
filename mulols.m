function beta = mulols(y, Y)
%*************************************************************************
% This function computes the multiple OLS estimator
%
%   INPUTS:
%       y : (n x p) data matrix, where
%            p is the number of time series and n is the length of them
%       Y : (n x nreg) matrix with regression variables, where
%            nreg is the number of regression variables
%
%   OUTPUT:
%    beta : (nreg x p) vector of regression coefficients
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
%*************************************************************************

[n, m] = size(y);
[mY, nY] = size(Y);
if n ~= mY
    error('y and Y must have the same number of rows in mulols');
end
%[Q,R]=qr(Y,0);
%qy=Q'*y;
%beta=R(1:nY,:)\qy(1:nY);
beta = Y \ y;
