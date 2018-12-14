function resid = var_res(y, nlag, x)
% PURPOSE: performs vector autogressive estimation
%          and returns only residuals
%---------------------------------------------------
% USAGE: resid = var_res(y,nlag,x)
% where:    y    = an (nobs x neqs) matrix of y-vectors
%           nlag = the lag length
%           x    = optional matrix of variables (nobs x nx)
%                 (NOTE: constant vector automatically included)
%---------------------------------------------------
% RETURNS: a matrix of residuals (nobs x neqs)
%---------------------------------------------------
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

[nobs, neqs] = size(y);

nx = 0;

if nargin == 3
    [nobs2, nx] = size(x);
    if (nobs2 ~= nobs)
        error('var_res: nobs in x-matrix not the same as y-matrix');
    end;
end;

% adjust nobs to feed the lags
nobse = nobs - nlag;

ylag = glags(y, nlag);

% form x-matrix
if nx
    xmat = [ylag, x(nlag+1:nobs, :), ones(nobse, 1)];
else
    xmat = [ylag, ones(nobse, 1)];
end;

beta = mulols(y(nlag+1:nobs, :), xmat);
resid = y(nlag+1:nobs, :) - xmat * beta;
