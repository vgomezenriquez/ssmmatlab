function xlag = ltflag(x, n)
% this function generates a matrix of n+1 lags from a matrix or vector to use
%          in ARMAX models. The first lag is the zero lag.
%---------------------------------------------------
% USAGE:     xlag = ltflag(x,nlag)
% where: x = a matrix of dimension nobs x nvar
%     nlag = number of lags to use in the ltf method
%---------------------------------------------------
% RETURNS:
%        xlag = a matrix of lags (nobs-nlag) x ((nlag+1)nvar)
%        x(t), x(t-1), x(t-2), ... x(t-nlag)
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

if nargin > 2
    error('ltflag: Wrong # of input arguments');
end;

[nobs, nvar] = size(x);

xlag = zeros(nobs-n, (n + 1)*nvar);
xlag(:, 1:nvar) = x(n+1:nobs, 1:nvar);
for i = 1:n
    xlag(:, 1+i*nvar:(i + 1)*nvar) = x(n-i+1:nobs-i, 1:nvar);
end
