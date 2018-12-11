function ycor = genycor(y, Y, ny, g)
%
% this function computes the series corrected by regression effects in a
% regression model of the form
%
% y = Y*g +e,
%
% where ny is the series length.
%
% Copyright (c) 21 July 2015 by Victor Gomez
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

[ng, mg] = size(g);
if ng > 0
    ycor = y - Y(1:ny, :) * g;
else
    ycor = y;
end
