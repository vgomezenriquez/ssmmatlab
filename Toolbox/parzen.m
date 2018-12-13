function [w, m] = parzen(n)
%
%        This function computes the weights for the
%        Parzen window
%
%     INPUTS:
%         n : lentgh of the series
%
%    OUTPUTS:
%         w : weights of the Parzen window
%         m : window lag size;
%
% Copyright (c) 21 July 2003 by Victor Gomez
% Ministerio de Hacienda, Direccion Gral. de Presupuestos,
% Subdireccion Gral. de Analisis y P.E.,
% Alberto Alcocer 2, 1-P, D-34, 28046, Madrid, SPAIN.
% Phone : +34-915835439
% E-mail: vgomez@sepg.minhap.es
%
% The author assumes no responsibility for errors or damage resulting from the use
% of this code. Usage of this code in applications and/or alterations of it should
% be referenced. This code may be redistributed if nothing has been added or
% removed and no money is charged. Positive or negative feedback would be appreciated.
%
m = floor(n/3-1);
w = zeros(m, 1);
for i = 1:floor(m/2)
    w(i) = 1 - (6 * (i^2) / m^2) * (1 - (i / m));
end
for i = floor(m/2) + 1:m
    w(i) = 2 * ((1 - (i / m))^3);
end
