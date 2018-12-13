function [gamser, gln] = gser(A, X)
%
% auxiliary function called in function gammp
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
gln = gammln(A);
if (X <= 0)
    %        if (X < 0) return; end
    gamser = 0;
    return
end
ap = A;
sum = 1 / A;
del = sum;
for n = 1:500
    ap = ap + 1;
    del = del * X / ap;
    sum = sum + del;
    if (abs(del) < abs(sum) * eps)
        break
    end
end
gamser = sum * exp(-X+A*log(X)-gln);
