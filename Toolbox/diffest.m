function [yd, beta] = diffest(y, Y, s, S, dr, ds, dS, est)
%*************************************************************************
% This function computes differences of the data vector and the matrix with
% regression variables and optionally performs estimation of regression
% coefficients
%
%  INPUTS:
%      y : data vector
%      Y : matrix with regression variables
%      s : frequency of the data, number of seasons
%      S : number of periods in each season
%     dr : regular differences
%     ds : seasonal differences corresponding to s (1 - B^s)^ds
%     dS : differences corresponding to S (1 - B^S)^dS
%    est = 1 : estimation of regression coefficients
%        = 0 : no estimation of regression coefficients
%
%  OUPUTS:
%     yd : data vector after differencing
%   beta : OLS estimator of the coefficients of the regression variables
%
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

[mY, nY] = size(Y);
n = length(y);
if nY > 0
    yd = [y, Y(1:n, :)];
else
    yd = y;
end
if dr > 0
    for i = 1:dr
        yd = diferm(yd, 1);
    end
end
if ds > 0
    yd = diferm(yd, s);
end
if dS > 0
    yd = diferm(yd, S);
end
if est == 1 & nY > 0
    beta = bols(yd(:, 1), yd(:, 2:nY+1));
else
    beta = [];
end
