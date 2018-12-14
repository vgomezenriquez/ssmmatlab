function [trendf, strendf, cicf, scicf] = stfilt4(Xt, Pt, xx, y)
%*************************************************************************
% This function copies the filtered states to the arrays
% trendf, cicf, and their Mse to the arrays strendf and scicf.
%
% Auxiliary function called in usa4vcv_d.m.
% Reference: ``Estimating Potential Output, Core Inflation
% and the NAIRU as Latent Variables'', by Rafael Domenech
% and Victor Gomez, Journal of Business and Economic Statistics (2006)
%
%    INPUTS:
%        Xt : an (n x nalpha matrix) containing the estimated x_{t|t}
%       Pt : an ((nalpha*n) x nalpha) matrix containing the
%            Mse of x_{t|t}
%       xx : vector with all parameters after estimation
%        y : data matrix
%
%   OUTPUTS:
%   trendf : filtered trend components of output, inflation, unemployment
%            and investment
%  strendf : standard errors of trendf
%     cicf : filtered cyclical components of output, inflation, unemployment
%            and investment
%    scicf : standard errors of cicf
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
trendf = zeros(nm, m);
cicf = zeros(nm, m);
strendf = zeros(nm, m);
scicf = zeros(nm, m);
for i = 1:nm
    ia = (i - 1) * mx + 1:i * mx;
    P = Pt(ia, 1:mx);
    xt = Xt(i, 1:mx);
    trendf(i, 1) = xt(1);
    strendf(i, 1) = sqrt(P(1, 1));
    trendf(i, 2) = xt(2);
    strendf(i, 2) = sqrt(P(2, 2));
    trendf(i, 3) = xt(3);
    strendf(i, 3) = sqrt(P(3, 3));
    trendf(i, 4) = xt(4);
    strendf(i, 4) = sqrt(P(4, 4));
    cicf(i, 1) = xt(mx);
    scicf(i, 1) = sqrt(P(mx, mx));
    cicf(i, 2) = xx(14) * xt(mx) + xx(15) * xt(mx-1) + xx(16) * xt(mx-2);
    scicf(i, 2) = sqrt(xx(16:-1:14)*P(mx-2:mx, mx-2:mx)*xx(16:-1:14)');
    cicf(i, 3) = xx(17) * xt(mx) + xx(18) * xt(mx-1);
    scicf(i, 3) = sqrt(xx(18:-1:17)*P(mx-1:mx, mx-1:mx)*xx(18:-1:17)');
    cicf(i, 4) = xx(19) * xt(mx);
    scicf(i, 4) = sqrt(P(mx, mx)*xx(19)^2);
end
