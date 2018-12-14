function [C, ierror] = pmatmul(A, B)
%
% This function computes the product of the polynomial matrices A and B
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
[na, ma, p] = size(A);
[nb, mb, q] = size(B);
C = [];
if (ma ~= nb)
    ierror = 1;
    return
end
iA = 0;
iB = 0;
if (na == ma) && ~any(any(A(:, :, 1)-eye(na)))
    iA = 1;
end
if (nb == mb) && ~any(any(B(:, :, 1)-eye(nb)))
    iB = 1;
end
C = zeros(na, mb, p+q-1);
for i = 1:p
    for j = 1:q
        if (i == 1) && (iA == 1)
            C(:, :, j) = C(:, :, j) + B(:, :, j);
        elseif (j == 1) && (iB == 1)
            C(:, :, i) = C(:, :, i) + A(:, :, i);
        else
            C(:, :, i+j-1) = C(:, :, i+j-1) + A(:, :, i) * B(:, :, j);
        end
    end
end
