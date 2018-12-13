function [gammcf, gln] = gacf(A, X)
%
% auxiliary function called in function gammp (gamma density function)
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
%


eps = 3.e-7;
fpmin = 1.e-30;
gln = gammln(A);
b = X + 1 - A;
c = 1 / fpmin;
d = 1 / b;
h = d;
%the following line added on 30-10-2008
gammcf = -1;
%if gammcf=-1 on output, the routine did not work
for i = 1:100
    an = -i * (i - A);
    b = b + 2;
    d = an * d + b;
    if (abs(d) < fpmin), d = fpmin;
    end
    c = b + an / c;
    if (abs(c) < fpmin), c = fpmin;
    end
    d = 1 / d;
    del = d * c;
    h = h * del;
    if (abs(del-1) < eps), gammcf = exp(-X+A*log(X)-gln) * h;
        break;
    end
end
