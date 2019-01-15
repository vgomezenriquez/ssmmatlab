function parm = tfivparm(Y, parm)
%
% this function adds to the structure parm the fields ninput and inputv
%
% Input arguments:
% Y    : a matrix containing the input variables
% parm : a structure where
%  s:  seasonality
%  S:  second seasonality
% .p:  AR order
% .ps: order of the AR of order s
% .q:  order of the regular MA
% .qs: order of the MA of order s (1 at most)
% .qS: order of the MA of order S (1 at most)
% .dr: order of regular differencing
% .ds: order of differencing of order s
% .dS: order of differencing of order S
% .pvar:  array containing the indices of variable parameters
% .pfix:  array containing the indices of fixed parameters
%
% Output arguments:
% parm: the input structure with the added fields
% .ninput: number of inputs
% .inputv: array containing the input variables
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

[nY, ninput] = size(Y);
parm.ninput = ninput;
parm.inputv = zeros(nY, ninput);
for i = 1:ninput
    parm.inputv(:, i) = Y(:, i);
end
