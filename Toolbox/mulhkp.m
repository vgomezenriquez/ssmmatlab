function [HKp] = mulhkp(H, K)
%
% this function multiplies matrix H by matrix K'
% assuming H*K' is symmetric
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

[r, n] = size(H);
HKp = size(r, r);
for i = 1:r - 1
    HKp(i:r, i) = H(i:r) * K(i);
    HKp(i, i+1:r) = HKp(i+1:r, i)';
end
HKp(r, r) = H(r) * K(r);
