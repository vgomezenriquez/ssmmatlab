function chk = chkroots(x, p, ps, q, qs, qS)
%*************************************************************************
% This function tests whether all roots of the polynomials of a
% multiplicative ARMA model are outside of the unit circle
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
%  OUTPUT:
%    chk = 0 : roots are outside the unit circle
%        = 1 : roots are not outside the unit circle
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

y = arpar(x, p, ps, q, qs, qS);
n = length(y);
chk = 0;
for i = 1:n
    if abs(y(i)) > 1
        chk = 1;
        break
    end
end
