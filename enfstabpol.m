function polt = enfstabpol(pol)
%
% This function multiplies the parameters of a matrix polynomial by some
% factor of the form .95^n until it becomes stable (all roots outside the
% unit circle). That is, Polynomial P(z) is transformed into P(lambda*z)
%
% Input arguments:
% pol: a matrix polynomial
%
% Output arguments:
% polt: where polt(z) = pol(lambda*z)
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

[n, m, k] = size(pol);
r = k - 1;
polt = pol;
if r == 0
    return
else
    th(:, :, 1) = eye(n);
    [H, F] = qarmax2ss2(pol, th);
    maeval = max(abs(eig(F)));
    epsilon = 0.015;
    c = 1. + epsilon;
    lambda = 1 / (c * maeval);
    li = lambda;
    for i = 1:r
        polt(:, :, i+1) = polt(:, :, i+1) * li;
        li = li * lambda;
    end
end
