function y = gammp(A, X)
%
% computes the probability density function of the gamma distribution
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


if (X < 0) | (A < 0)
    if (A < 0), y = 2;
        return;
    end
end
if (X < A + 1)
    [gamser, gln] = gser(A, X);
    y = gamser;
else
    [gammcf, gln] = gacf(A, X);
    y = 1 - gammcf;
end
