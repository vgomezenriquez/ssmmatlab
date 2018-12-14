function X = glags(x, lags)
%
% this function generates a matrix containing lags of x
%
% Input arguments:
% x: an (n) x (1) vector series
% lags: the number of lags to be generated
%
% Output argument:
% X: an (n-lags) x (lags) matrix containing the lagged variables
% Note that lags observations are lost
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

[n, junk] = size(x);
X = [];
for i = 1:lags
    X = [X, x(lags-i+1:n-i, :)];
end
