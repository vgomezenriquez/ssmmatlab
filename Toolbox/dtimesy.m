function [yd, ferror] = dtimesy(D, y)
%
%
% This function obtains the series
%         yd_t = D(B)*y_t,
% where D(z)=D_0 + D_1*z + .... + D_r*z^r is a polynomial matrix compatible
% with y_t and B is the backshift operator, By_t = y_{t-1}.
%
% Input arguments:
%                  D: a n x m x r polynomial matrix
%                  y: an m x s matrix
% Output arguments:
%                  yd: the series D(B)*y_t
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
yd = [];
ferror = 0;
[nd, md, rd] = size(D);
[n, m] = size(y);
r = rd - 1;

if (m ~= md)
    disp('dimension mismatch in dtimesy')
    ferror = 1;
    return
end

yd = y(r+1:n, :) * D(:, :, 1)';
for i = 1:r
    yd = yd + y(r-i+1:n-i, :) * D(:, :, i+1)';
end
