function x = inest(yd, beta, s, S, p, ps, q, qs, qS, ols, a)
%
% this function computes initial ARMA estimates using
% the Hannan-Rissanen method
%
% Input arguments:
% yd   : an (n x m) matrix containing the series, yd(:,1), and an (n x m-1)
%        matrix of regression variables if m > 1.
% beta : an m-1 vector containing the OLS estimators if m > 1, empty if m =
%        1
% s    : seasonality
% S    : second seasonality
% p    : degree of AR polynomial
% ps   : degree of AR seasonal polynomial
% q    : degree of MA polynomial
% qs   : degree of MA seasonal polynomial
% qS   : degree of MA second seasonal polynomial
% ols  : = 1, perform OLS, = 0, use the Durbin Levinson algorithm in the HR
%        method
% a    : an integer, the degree of the AR approximation in the first step
%        of the Hanna-Rissanen method.
%
% Output arguments:
% x    : a vector containing the ARMA parameter estimates
%
%
% Copyright (c) 21 July 2014 by Victor Gomez
% Ministerio de Hacienda y A.P., Direccion Gral. de Presupuestos,
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

[n, m] = size(yd);
% correct series for regression effects
if m > 1
    yc = yd(:, 1) - yd(:, 2:m) * beta;
else
    yc = yd(:, 1);
end
x = hanris(yc, s, S, p, ps, q, qs, qS, ols, a);
% filter the series by the autoregressive polynomial
% to obtain better MA estimates
if (p + ps > 0) && (q + qs + qS > 0)
    yfp = zeros(n-p, 1);
    if p > 0
        for i = p + 1:n
            yfp(i-p) = [1, x(1:p)] * yc(i:-1:i-p);
        end
    elseif p == 0
        yfp = yc;
    end
    yfpps = zeros(n-p-ps*s, 1);
    if ps > 0
        for i = ps * s + 1:n - p
            yfpps(i-ps*s) = [1, x(p+1:p+ps)] * yfp(i:-s:i-ps*s);
        end
    elseif ps == 0
        yfpps = yfp;
    end
    xx = hanris(yfpps, s, S, 0, 0, q, qs, qS, ols, a);
    nx = length(x);
    x(p+ps+1:nx) = xx;
end
