function jpr = jprod(x, y, p, q)
%      given two n-column vectors, x and y, and a signature matrix J=diag(I_p,-I_q),
%      this function calculates the product x'*J*y.
%
%---------------------------------------------------
% USAGE: jpr = jprod(x,y,p,q)
% where:    x, y = n-column vectors
%           p,q= integers such that J = diag(I_p,-I-q) is a signature
%           matrix
%---------------------------------------------------
% RETURNS: the J product of x and y
%---------------------------------------------------
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
%
[n, junk] = size(x);
[m, junk] = size(y);
if n ~= m
    error('x and y have not the same length in jprod')
end
if n ~= p + q
    error('n nonequal to p + q in jprod')
end
jpr = 0.d0;
if p > 0
    jpr = jpr + x(1:p)' * y(1:p);
end
if q > 0
    jpr = jpr - x(p+1:p+q)' * y(p+1:p+q);
end
