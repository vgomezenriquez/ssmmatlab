function y = parar(x, p, ps, q, qs, qS)
%*************************************************************************
% Given the polynomials of a multiplicative ARMA model such that their
% coefficients have been transformed into parcor coefficients, this
% function converts parcor coefficients back into AR coefficients
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
%      y : coefficients of all the polynomials of the
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

y = zeros(size(x));

if p > 0
    y(1:p) = pafi(x(1:p));
end
if ps > 0
    y(p+1:p+ps) = pafi(x(p+1:p+ps));
end
if q > 0
    y(p+ps+1:p+ps+q) = pafi(x(p+ps+1:p+ps+q));
end
if qs > 0
    y(p+ps+q+1:p+ps+q+qs) = pafi(x(p+ps+q+1:p+ps+q+qs));
end
if qS > 0
    y(p+ps+q+qs+1:p+ps+q+qs+qS) = pafi(x(p+ps+q+qs+1:p+ps+q+qs+qS));
end
