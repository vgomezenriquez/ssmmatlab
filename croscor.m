function [cr, stdx, stdy] = croscor(x, y, lag)
%
%
%        This function computes the correlation between x(t) and
%        y(t+lag). Lag can be both positive and negative.
%
%     INPUTS:
%------------
%    REQUIRED
%       x,y : series; y = x, if autocorrelations are to be computed
%
%    OPTIONAL
%       lag : number of lags at which the correlations are to be computed
%             lag = length(y)-1, if lag is not input to croscor or if lag
%             is empty
%
%    OUTPUTS:
%------------
%        cv : correlations between x and y
%      stdx : standard deviation of x
%      stdy : standard deviation of y
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
n = length(y);
if n ~= length(x)
    error('x and y must be the same size.');
end
stdx = sqrt((x - mean(x))'*(x - mean(x)));
stdy = sqrt((y - mean(y))'*(y - mean(y)));
if lag < 0
    cr = ((x(1-lag:n) - mean(x))' * (y(1:n+lag) - mean(y))) / (stdx * stdy);
elseif lag > 0
    cr = ((x(1:n-lag) - mean(x))' * (y(1+lag:n) - mean(y))) / (stdx * stdy);
else
    cr = ((x - mean(x))' * (y - mean(y))) / (stdx * stdy);
end
