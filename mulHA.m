function [C, ierror] = mulHA(H, A, kro)
%
% This function computes the product of the polynomial matrices H and A
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
[nh, mh] = size(H);
C = [];
if (mh ~= na)
    ierror = 1;
    return
end
C = zeros(nh, ma);
r = 0;
for i = 1:nh
    if kro(i) == 0
        C(i, :) = H(i, :) * A;
    else
        C(i, :) = A(r+1, :);
    end
    r = r + kro(i);
end
