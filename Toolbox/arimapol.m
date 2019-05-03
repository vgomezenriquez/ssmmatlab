function [phirs, alprsS, thrsS] = arimapol(x, s, S, p, ps, dr, ds, dS, q, qs, qS)
%
% this function computes the different polynomials for an ARIMA model
%
%     INPUTS:
%           x:  an array containing the ARIMA parameter values
%           s:  seasonality
%           S:  second seasonality
%           p:  AR order
%          ps:  order of the AR of order s
%           q:  order of the regular MA
%          qs:  order of the MA of order s (1 at most)
%          qS:  order of the MA of order S (1 at most)
%          dr: order of regular differencing
%          ds: order of differencing of order s
%          dS: order of differencing of order S
%  OUTPUTS:
%     phirs  : an array containing the AR polynomial
%     alprsS : an array containing the differencing polynomial
%     thrsS  : an array containing the MA polynomial
%
% Copyright (c) 21 July 2015 by Victor Gomez
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
alpr = 1;
if dr > 0
    delta = [-1, 1];
    for i = 1:dr
        alpr = conv(delta, alpr);
    end
end
alps = 1;
if ds > 0
    delta = [-1, 1];
    K = ones(1, s);
    alps = conv(delta, K);
end
alpS = 1;
if dS > 0
    delta = [-1, 1];
    K = ones(1, S);
    alpS = conv(delta, K);
end
alprs = conv(alpr, alps);
alprsS = conv(alprs, alpS);
%
phi = 1;
phis = 1;
th = 1;
ths = 1;
thS = 1;
if p > 0
    phi = conv(phi, [fliplr(x(1:p)), 1]);
end
if ps > 0
    phis = zeros(1, s * ps +1);
    for i = 1:ps
       phis(1 + s * (i - 1)) = x(p + ps + 1 - i);
    end
    phis(s * ps + 1) = 1;
end
phirs = conv(phi, phis);
if q > 0
    th = conv(th, [fliplr(x(p+ps+1:p+ps+q)), 1]);
end
if qs > 0
    ths = zeros(1, s * ps +1);
    for i = 1:qs
       ths(1 + s * (i - 1)) = x(p + ps + q + qs + 1 - i);
    end
    ths(s * qs + 1) = 1;
end
if qS > 0
    thS = zeros(1, S+1);
    thS(1) = x(p+ps+q+qs+qS);
    thS(S+1) = 1;
end
thrs = conv(th, ths);
thrsS = conv(thrs, thS);
