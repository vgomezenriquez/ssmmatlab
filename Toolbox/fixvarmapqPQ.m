function [str, ferror] = fixvarmapqPQ(str)
% PURPOSE: given a structure containing information about a VARMA model, it
% fixes the parameters in phi, Phi, th, Th and Lh that correspond to those
% elements in phin, Phin, thn, Thn and Lhn that are not equal to NaN.
%---------------------------------------------------
% USAGE: str = fixvarmaxpqPQ(str)
% where:   str      = a structure created with function suvarmapqPQ
%---------------------------------------------------
% RETURNS: str = a structure containing model information
%---------------------------------------------------
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

ferror = 0;

if (~isfield(str, 'phin') | ~isfield(str, 'Phin') | ~isfield(str, 'thn') | ...
        ~isfield(str, 'Thn') | ~isfield(str, 'Lhn'))
    ferror = 1;
    disp('some of the fields phin, Phin, thn, Thn or Lhn are not in fixvarmapqPQ');
    return
end

phin = str.phin;
thn = str.thn;
Phin = str.Phin;
Thn = str.Thn;
Lhn = str.Lhn;

[np, mp, pm1] = size(phin);
[nt, mt, qm1] = size(thn);
[nP, mP, Pm1] = size(Phin);
[nT, mT, Qm1] = size(Thn);

%create vector for fixed model parameters
x = [];
for i = 2:pm1
    x = [x, vec(phin(:, :, i))'];
end
for i = 2:Pm1
    x = [x, vec(Phin(:, :, i))'];
end
for i = 2:qm1
    x = [x, vec(thn(:, :, i))'];
end
for i = 2:Qm1
    x = [x, vec(Thn(:, :, i))'];
end
x = [x, Lhn(2:end)'];

nparm = str.nparm;
xi = str.xi;
xx = str.x;
xv = [];
xf = [];
for i = 1:nparm
    if isnan(x(i))
        xv = [xv, xx(i)];
    else
        xf = [xf, xx(i)];
        xi(i) = 0;
    end
end
str.xi = xi;
str.xv = xv;
str.xf = xf;