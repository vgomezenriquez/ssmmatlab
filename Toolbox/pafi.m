function fi = pafi(r)
%
% this function transforms the partial autocorrelation coefficients
% of an autoregressive model 1+phi_1*z+...+phi_p*z^p into the
% model parameters
% input : r, a 1 x p vector
% output: fi, a 1 x p vector
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

p = length(r);
fi = r;
if p == 1
    return
else
    for i = 2:p
        fi(1:i-1) = fi(1:i-1) + fi(i) * fliplr(fi(1:i-1));
    end
end
