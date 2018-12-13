function [b, a, err] = shank(nu, db, da)
%
% Given the impulse response function nu(z) = nu_0 + nu_1*z + nu_2*z^2 +
% ...., this function obtains the rational approximation nu(z ) ~=
% b(z)/a(z), where b(z) = b_0 + b_1*z + ... + b(nb)*z^nb and a(z) = 1
% + a_1*z + ... + a(na)*z^na. The coefficients of b(z) and a(z) are
% computed using Shank's method (Discrete Random signal processing, p. 558)
%
% Input arguments:
% nu: the impulse response function
% db: the degree of b(z)
% da: the degree of a(z)
%
% Output arguments:
% b: a (db+1) array containing the b(z) coefficients in ascending order
% a: a (da+1) array containing the a(z) coefficients in ascending order
%
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

nnu = length(nu);
err = 0;
r = max(db, da);
if (nnu - r - 1 < da)
    disp('not enough nu weights in shank')
    err = 1;
    return
end
a = ones(da+1, 1);
b = zeros(db+1, 1);
Ha = [];
h = zeros(1, da+1);
for i = r + 2:nnu
    for j = 1:da + 1
        h(j) = nu(i-j+1);
    end
    Ha = [Ha; h];
end
% Ha
if isempty(Ha)
    b = nu; %da is zero and nnu = db + 1
    %  a,b,pause
    return
end
a(2:end) = Ha(:, 2:end) \ (-Ha(:, 1));
ha = poldiv(1, a, nnu-1);
HA = zeros(nnu, db+1);
hb = zeros(nnu, 1);
for j = 1:nnu, hb(j) = nu(j);
end
for j = 1:db + 1
    HA(j:end, j) = ha(1:end-j+1);
end
% HA
% hb
b = HA \ hb;
