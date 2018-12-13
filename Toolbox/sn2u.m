function [p, ierrsn2u] = sn2u(c)
% This function transforms a polynomial in the variables S(n)=z^n + z^(-n)
% into a polynomial in the variable U=z + z^(-1).
% The recursions
%   S(n) = U*S(n-1) - S(n-2),    n>=2
%   S(0)=2,  S(1)=U,
% are used.
%
% c on input contains the coefficients in the order 1,2,...,n,
%  c= c(1) + c(2)*S(1) + ... + c(n-1)*S(n-2) + c(n)*S(n-1)
% p on output contains the coefficients in reversed order n, n-1, ... ,1,
%  p(U) = p(1)*U^(n-1) + p(2)*U^(n-2) + ... + p(n-1)*U + p(n)
%
%
% Copyright (c) 21 July 2003 by Victor Gomez
% Ministerio de Hacienda, Direccion Gral. de Presupuestos,
% Subdireccion Gral. de Analisis y P.E.,
% Alberto Alcocer 2, 1-P, D-34, 28046, Madrid, SPAIN.
% Phone : +34-915835439
% E-mail: vgomez@sepg.minhac.es
%
% The author assumes no responsibility for errors or damage resulting from the use
% of this code. Usage of this code in applications and/or alterations of it should
% be referenced. This code may be redistributed if nothing has been added or
% removed and no money is charged. Positive or negative feedback would be appreciated.


n = length(c);
p = [];
ierrsn2u = 0;
if (n == 1)
    p = c(1);
elseif (n >= 2)
    p = [c(2), c(1)];
    Sn1 = [1, 0];
    U = [1, 0];
    Sn2 = [0, 0, 2];
    for i = 3:n
        S = conv(Sn1, U) - Sn2;
        q = c(i) * S;
        p0 = [0, p];
        p = p0 + q;
        Sn2 = [0, 0, Sn1];
        Sn1 = S;
    end
end