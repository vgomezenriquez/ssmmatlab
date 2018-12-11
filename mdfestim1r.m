function [D, DA, yd, ferror] = mdfestim1r(y, x, prt, nr)
%
% This function performs VARMAX(1,1) estimation with rank imposed and
% returns 'differencing' polynomial and 'differenced' series.
%
% Inputs: y: matrix containing the output series
%         x: matrix containing the input series
%         prt      = 1 print results of the VARX, VARMAX(p,p,p) tests
%         nr: an integer, number of regular unit roots
%  Output: D: an (ny x ny x 2) 'differencing' matrix polynomial
%         yd= matrix containing the 'differenced' series
%         DA= matrix of the form [DAr Indxr], where DAr is the
%             parameterization of the differencing polynomial, and Indxr is
%             an index vector to identify the l.i. rows of DAr.
%         ferror= a flag for errors
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

yd = [];
D = [];
DA = [];

[ny, s] = size(y);
if (s == 1)
    ferror = 2;
    disp('mdfestim1: number of series should be greater than one')
    return
end
[nx, m] = size(x);
if ~isempty(x)
    if (nx ~= ny)
        ferror = 3;
        disp('mdfestim1: nobs in x-matrix not the same as y-matrix');
        return
    end;
end

seas = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%determination of p
hr3 = 1; %perform only the first two stages of the HR method
maxorder = 0;
ct = [];
[order0, kro0] = varmaxgenid(y, x, seas, maxorder, hr3, ct, prt);
lagsopt = kro0(1);
if prt == 1
    out = lagsopt;
    fprintf(1, 'Estimated p in VARMAX(p,p,p): ')
    fprintf(1, 'p = %2d\n', out);
end
p = lagsopt;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initial estimation using all Kronecker indices equal to p
kro = repmat(p, 1, s);
str = estvarmaxkro(y, x, seas, kro, hr3);
residv = str.residv; %first residuals obtained with VARMAX(p,p,p)
resid = str.resid2;
nres = size(resid, 1);
if ny > nres
    residv = [residv(1:ny-nres, :); resid];
else
    residv = resid;
end
a = residv;
clear str
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%estimate the model with the rank r imposed
r = s - nr; %rank of Pi1
ir = 1;
[Pi1r, alpha, betap, betaor, ferror] = evarma11rurimp(y, x, a, ir, r);
% 'differencing' polynomial
D(:, :, 1) = eye(s);
D(:, :, 2) = -betaor * pinv(betaor'*betaor) * betaor';
% 'differenced' series
[yd, ferror] = dtimesy(D, y);
%parameterization
[DA, ferror] = parambeta(betaor);
