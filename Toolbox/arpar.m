function y = arpar(x, p, ps, q, qs, qS)
%*************************************************************************
%
% Given the polynomials of a multiplicative ARMA model, this function
% transforms the coefficients of each polynomial as if they were AR
% coefficients into partial correlation coefficients. This may be used as a
% test for stationarity because a polynomial is stable if, and only if, all
% its partial coefficients are less than one in absolute value.
%
%  INPUTS:
%      x : coefficients of the polynomials of a multiplicative ARMA model
%      p, ps, q, qs, qS : integers specifying where the coefficients of the
%      ARMA model are in x.
%      More specifically,
%      p : first p are AR coefficients
%     ps : starting with the (p+1)th coefficient, the next ps are AR c.
%      q : starting with the (p+1+ps+1)th coefficient, the next q are MA c.
%     qs : starting with the (p+1+ps+1+q+1)th coefficient,
%          the next qs are MA c.
%     qS : starting with the (p+1+ps+1+q+1+qs+1)th coefficient,
%          the next qS are MA c.
%
% OUTPUTS:
%      y : partial correlation coefficients of all the polynomials of the
%          ARMA model
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

pps = p + ps;
ppsq = pps + q;
ppsqqs = ppsq + qs;
ppsqqsqS = ppsqqs + qS;
y = zeros(size(x(1:ppsqqsqS)));

if p > 0
    y(1:p) = fipa(x(1:p));
end
if ps > 0
    y(p+1:pps) = fipa(x(p+1:pps));
end
if q > 0
    y(pps+1:ppsq) = fipa(x(pps+1:ppsq));
end
if qs > 0
    y(ppsq+1:ppsqqs) = fipa(x(ppsq+1:ppsqqs));
end
if qS > 0
    y(ppsqqs+1:ppsqqsqS) = fipa(x(ppsqqs+1:ppsqqsqS));
end
