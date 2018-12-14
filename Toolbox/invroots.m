function z = invroots(x, p, ps, q, qs, qS)
%*************************************************************************
% This function enforces that all roots of the polynomials of a
% multiplicative ARMA model are outside of the unit circle. This is
% achieved by first transforming the coefficients of each polynomial as if
% they were AR coefficients into partial correlation coefficients. Then, if
% there are parcor coefficients greater than one in absolute value, these
% are transformed into parcor coefficients that lie in (-1,1). Finally, the
% new parcor coefficients are transformed into AR coefficients.
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
%      z : coefficients such that all roots are outside of the unit circle
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


z = arpar(x, p, ps, q, qs, qS);
n = length(z);
for i = 1:n
    zi = z(i);
    if (zi < -1.) || (zi > 1.)
        %   z(i)=(exp(zi)-1.)/(exp(zi)+1.);
        z(i) = atan(zi) * 2 ./ pi;
    end
end
z = parar(z, p, ps, q, qs, qS);
