function [c0, cv, r] = autcov(y, lag, ic)
%*************************************************************************
%
% This function computes the sample autocovariances and autocorrelations of
% y(t) up to specified lag. The variable has been previously demeaned.
%
%   INPUTS:
%       y : input vector
%     lag : integer specifying up to which lag cv and/or r are computed
%      ic = 1: compute autocorrelations
%           0: do not compute autocorrelations
%
%  OUTPUTS:
%      c0 : variance of y
%      cv : autocovariances of y
%       r : autocorrelations of y
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

n = length(y);
cv = zeros(lag, 1);
r = [];
% variance
c0 = y' * y / n;
% autocovariances
for i = 1:lag
    cv(i) = y(1:n-i)' * y(1+i:n) / n;
end
% autocorrelations
if ic == 1
    r = cv / c0;
end
