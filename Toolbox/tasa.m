function yt = tasa(y, s)
% ************************************************************************
%  This function computes growth rate of order s
%
%   INPUTS:
%       y : data vector
%       s : integer, order of the growth rate of y
%
%   OUTPUT:
%      yt : growth rate of y of order s
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
%*************************************************************************

n = length(y);
yt = (y(s+1:n) - y(1:n-s)) ./ y(1:n-s);
