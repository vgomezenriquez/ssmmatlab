function nr = mnorm(x)
%      given an n-vector x, this function calculates the square of the
%      euclidean norm of x.
%
%      the euclidean norm is computed by accumulating the sum of
%      squares in three different sums. the sums of squares for the
%      small and large components are scaled so that no overflows
%      occur. non-destructive underflows are permitted. underflows
%      and overflows do not occur in the computation of the unscaled
%      sum of squares for the intermediate components.
%      the definitions of small, intermediate and large components
%      depend on two constants, rdwarf and rgiant. the main
%      restrictions on these constants are that rdwarf^2 not
%      underflow and rgiant^2 not overflow. the constants
%      given here are suitable for every known computer.
%---------------------------------------------------
% USAGE: nr = mnorm(x)
% where:    x = an n-dimensional vector
%---------------------------------------------------
% RETURNS: the square of the euclidean norm of x
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
one = 1.0d0;
zero = 0.0d0;
rdwarf = 3.834d-20;
rgiant = 1.304d19;
s1 = zero;
s2 = zero;
s3 = zero;
x1max = zero;
x3max = zero;
floatn = double(n);
agiant = rgiant / floatn;
for i = 1:n
    xabs = abs(x(i));
    if ((xabs > rdwarf) & (xabs < agiant))
        %
        %    sum for intermediate components.
        %
        s2 = s2 + xabs^2;
    elseif (xabs <= rdwarf)
        %
        %   sum for small components.
        %
        if (xabs <= x3max)
            if (xabs ~= zero)
                s3 = s3 + (xabs / x3max)^2;
            end
        else
            s3 = one + s3 * (x3max / xabs)^2;
            x3max = xabs;
        end
    else
        %
        %   sum for large components.
        %
        if (xabs <= x1max)
            s1 = s1 + (xabs / x1max)^2;
        else
            s1 = one + s1 * (x1max / xabs)^2;
            x1max = xabs;
        end
    end
end
%
%  calculation of the square of the norm.
%
if (s1 == zero)
    if (s2 == zero)
        nr = (x3max^2) * s3;
    else
        if (s2 >= x3max)
            nr = s2 * (one + (x3max / s2) * (x3max * s3));
        else
            nr = x3max * ((s2 / x3max) + (x3max * s3));
        end
    end
else
    nr = (x1max^2) * (s1 + (s2 / x1max) / x1max);
end
