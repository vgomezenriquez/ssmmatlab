function [Morp, ferror] = m2mor(M)
%
% Given an (s x r+1) matrix M, this function obtains an (s-r x s) matrix
% Morp such that Morp*M = 0. The last column of M is assumed to be an index
% Idx, such that for the i-th row of M, Idx(i) = 0 means that this row is
% linearly independent.
%
% Inputs  :    M: an (s x r+1) matrix
%  Output : Morp: an (s-r x s) matrix such that Morp*M = 0.
%
% Copyright (c) 21 July 2014 by Victor Gomez
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

Morp = [];
ferror = 0;

[s, rp1] = size(M);
Idx = M(:, rp1);
r = rp1 - 1;
if sum(Idx) ~= s - r
    ferror = 1;
    return
end
Morp = zeros(s-r, s);
B = zeros(s-r, r);
Ilix = zeros(1, r);
Ildx = zeros(1, s-r);
cont = 0;
contli = 0;
for i = 1:s
    if (Idx(i) == 1)
        cont = cont + 1;
        B(cont, :) = -M(i, 1:r);
        Ildx(cont) = i;
    else
        contli = contli + 1;
        Ilix(contli) = i;
    end
end
Morp(:, Ildx) = eye(s-r);
Morp(:, Ilix) = B;
