function [trends, strends, cics, scics] = stsmt4(Xt, Pt, xx, y, sconp)
%*************************************************************************
% This function copies the smoothed states to the arrays
% trends, cics, and their Mse to the arrays strends and scics.
%
% Auxiliary function called in usa4vcv_d.m.
% Reference: ``Estimating Potential Output, Core Inflation
% and the NAIRU as Latent Variables'', by Rafael Domenech
% and Victor Gomez, Journal of Business and Economic Statistics (2006)
%
%    INPUTS:
%       Xt : an (n x nalpha matrix) containing the estimated x_{t|n}
%       Pt : an ((nalpha*n) x nalpha) matrix containing the
%            Mse of x_{t|n}
%       xx : vector with all parameters after estimation
%        y : data matrix
%    sconp : square root of concentrated parameter
%
%   OUTPUTS:
%   trends : smoothed trend components of output, inflation, unemployment
%            and investment
%  strends : standard errors of trends
%     cics : smoothed cyclical components of output, inflation, unemployment
%            and investment
%    scics : standard errors of cics
%
% Copyright (c) 21 July 2003 by Victor Gomez
% Ministerio de Hacienda, Direccion Gral. de Presupuestos,
% Subdireccion Gral. de Analisis y P.E.,
% Alberto Alcocer 2, 1-P, D-34, 28046, Madrid, SPAIN.
% E-mail: VGomez@sepg.minhap.es
%
% The author assumes no responsibility for errors or damage resulting from the use
% of this code. Usage of this code in applications and/or alterations of it should
% be referenced. This code may be redistributed if nothing has been added or
% removed and no money is charged. Positive or negative feedback would be appreciated.
%*************************************************************************


[n, m] = size(y);
[nx, mx] = size(Xt); %mx is the length of the state vector
nm = n - 4; %nm is the effective number of observations
trends = zeros(nm, m);
cics = zeros(nm, m);
strends = zeros(nm, m);
scics = zeros(nm, m);
for i = 1:nm
    ia = (i - 1) * mx + 1:i * mx;
    P = Pt(ia, 1:mx);
    xt = Xt(i, 1:mx);
    trends(i, 1) = xt(1);
    strends(i, 1) = sconp * sqrt(P(1, 1));
    trends(i, 2) = xt(2);
    strends(i, 2) = sconp * sqrt(P(2, 2));
    trends(i, 3) = xt(3);
    strends(i, 3) = sconp * sqrt(P(3, 3));
    trends(i, 4) = xt(4);
    strends(i, 4) = sconp * sqrt(P(4, 4));
    cics(i, 1) = xt(mx);
    scics(i, 1) = sconp * sqrt(P(mx, mx));
    cics(i, 2) = xx(14) * xt(mx) + xx(15) * xt(mx-1) + xx(16) * xt(mx-2);
    scics(i, 2) = sconp * sqrt(xx(16:-1:14)*P(mx-2:mx, mx-2:mx)*xx(16:-1:14)');
    cics(i, 3) = xx(17) * xt(mx) + xx(18) * xt(mx-1);
    scics(i, 3) = sconp * sqrt(xx(18:-1:17)*P(mx-1:mx, mx-1:mx)*xx(18:-1:17)');
    cics(i, 4) = xx(19) * xt(mx);
    scics(i, 4) = sconp * sqrt(P(mx, mx)*xx(19)^2);
end
