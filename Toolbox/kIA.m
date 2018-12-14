function [C, ierror] = kIA(A)
%
% This function computes the Kronecker product of the matrix I by A
%
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
%

ierror = 0;
[na, ma] = size(A);
k = na;
k2 = k * k;
C = zeros(k2);
if (na ~= ma)
    ierror = 1;
    return
end
for i = 1:k
    C((i - 1)*k+1:k*i, (i - 1)*k+1:k*i) = A;
end
