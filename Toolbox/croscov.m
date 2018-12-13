function cv = croscov(x, y, lag)
%
%    This function computes the covariance between x(t) and y(t+lag).
%    Lag can be both positive and negative.
%
%
%     INPUTS:
%------------
%    REQUIRED
%       x,y : series; y = x, if autocovariances are to be computed
%
%    OPTIONAL
%       lag : number of lags at which covariances are to be computed
%             lag = length(y)-1, if lag is not input to croscov or if lag
%             is empty
%
%     OUTPUT:
%------------
%        cv : covariances between x and y
%
% Copyright (c) 21 July 2003 by Victor Gomez
% Ministerio de Hacienda, Direccion Gral. de Presupuestos,
% Subdireccion Gral. de Analisis y P.E.,
% Alberto Alcocer 2, 1-P, D-34, 28046, Madrid, SPAIN.
% Phone : +34-915835439
% E-mail: vgomez@sepg.minhap.es
%
% The author assumes no responsibility for errors or damage resulting from the use
% of this code. Usage of this code in applications and/or alterations of it should
% be referenced. This code may be redistributed if nothing has been added or
% removed and no money is charged. Positive or negative feedback would be appreciated.
%
if nargin < 2 || isempty(x) || isempty(y)
    error('x and y are required inputs to croscov');
end
n = length(y);
if n ~= length(x)
    error('x and y must be the same size');
end
if nargin == 2 || isempty(lag)
    lag = n - 1;
end
if lag < 0
    cv = ((x((1 - lag):n) - mean(x))' * (y(1:(n + lag)) - mean(y))) / n;
elseif lag > 0
    cv = ((x(1:(n - lag)) - mean(x))' * (y((1 + lag):n) - mean(y))) / n;
else
    cv = ((x - mean(x))' * (y - mean(y))) / n;
end
