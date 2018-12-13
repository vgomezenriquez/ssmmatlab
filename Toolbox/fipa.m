function r = fipa(fi)
%
% this function transforms the parameters of an autoregressive
% model 1+phi_1*z+...+phi_p*z^p into its partial autocorrelation
% coefficients
% input : fi, a 1 x p vector
% output: r, a 1 x p vector
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

p = length(fi);
r = fi;
if p == 1
    return
else
    for i = p - 1:-1:1
        j = floor(i/2);
        s = fliplr(r(1:i));
        b = r(i+1);
        if mod(i, 2) == 1 & j > 0
            r(1:j) = (r(1:j) - b * s(1:j)) / (1 - b^2);
            r(j+1) = r(j+1) / (1 + b);
            r(j+2:i) = (r(j+2:i) - b * s(j+2:i)) / (1 - b^2);
        else
            r(1:i) = (r(1:i) - b * s) / (1 - b^2);
        end
    end
end
