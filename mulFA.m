function [C, ierror] = mulFA(F, A, kro)
%
% This function computes the product of the matrices F and A, where F is
% assumed to be in echelon form (x_{t+1} = F*x_{t} + K*a_{t}).
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
[nf, mf] = size(F);
C = [];
if (mf ~= na)
    ierror = 1;
    return
end
C = zeros(nf, ma);
j = 1;
%the following three lines added on 19-1-2011
while kro(j) == 0
    j = j + 1;
end
r = kro(j);
for i = 1:nf - 1
    if i == r
        C(i, :) = F(i, :) * A;
        j = j + 1;
        %the following three lines added on 19-1-2011
        while kro(j) == 0
            j = j + 1;
        end
        r = r + kro(j);
    else
        C(i, :) = A(i+1, :);
    end
end
C(nf, :) = F(nf, :) * A;
