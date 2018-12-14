function At = pmattrans(A)
%
% This function computes the transpose of a polynomial matrix A
%---------------------------------------------------
% USAGE: At = ptransmat (A)
% where:    A = a polynomial matrix not necessarily square
%---------------------------------------------------
% RETURNS:
%           At = the transpose of A
%---------------------------------------------------
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

K = [];
[s1, s2, np] = size(A);
At = zeros([s2, s1, np]);
for i = 1:np
    At(:, :, i) = A(:, :, i)';
end