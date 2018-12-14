function [beta, tv, d, M] = btval(x, ydf)
%*************************************************************************
%  This function computes the OLS estimator, the t-values and the standard
%  errors
%
%   INPUTS:
%       x : vector with parameters
%     ydf : matrix with the data and regression variables
%
%  OUTPUTS:
%    beta : OLS estimator
%      tv : t-value of beta
%       d : standard error of beta
%       M : covariance matrix of beta
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

[nd, md] = size(ydf);
md1 = md - 1;
[beta, M, e] = bmols(ydf(:, 1), ydf(:, 2:md));
vin = e' * e / (nd - md1 - length(x));
M = M * vin;
d = sqrt(diag(M));
% t-values
tv = beta ./ d;
