function [Lambda, alpha, betap, th, Th, L, ferror] = pr2ecf(xv, xf, DA, str)
% PURPOSE: given a structure containing information about a VARMA model in
% error correction form, it obtains the model
%---------------------------------------------------
% USAGE: [Lambda,alpha,betap,th,Th,L,ferror] = pr2ecf(xv,xf,DA,str)
% where:
%          xv       = a vector containing the parameters to be estimated
%          xf       = a vector containing the fixed parameters
%                     be estimated, =0, not
%          DA       = the matrix containing the parameterization of the
%                     differencing polynomial
%          str      = a structure containing the initial model information
%---------------------------------------------------
%---------------------------------------------------
% RETURNS:
%          Lambda   = a polynomial matrix of degree p-1, where p is the
%                     degree of the overall AR matrix polynomial
%          alpha    = an (s x r) matrix, where s is series dimension and r
%                     is the cointegration rank
%          betap    = an (r x s) matrix
%          th       = the regular MA matrix polynomial
%          Phi      = the seasonal AR matrix polynomial
%          Th       = the seasonal MA matrix polynomial
%          L        = the Cholesky factor of the innovations covariance
%                     matrix
%          ferror   = flag for errors
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
Lambda = [];
alpha = [];
betap = [];
th = [];
Th = [];
L = [];
if (~isfield(str, 'Lambda') || ~isfield(str, 'alpha') || ...
        ~isfield(str, 'betap') || ~isfield(str, 'th') || ~isfield(str, 'Th') || ...
        ~isfield(str, 'Lh'))
    ferror = 1;
    disp('some of the fields Lambda, alpha, betap, th, Th or Lh')
    disp('are not in pr2ecf');
    return
end
if (~isfield(str, 'xi'))
    ferror = 1;
    disp('the field xi is not in pr2ecf');
    return
end

xi = str.xi;
if isempty(xv)
    x = xf;
else
    if isempty(xf)
        x = xv;
    else
        x = zeros(size(xi));
        indvar = find(xi);
        indfix = find(~xi);
        x(indvar) = xv;
        x(indfix) = xf;
    end
end

%pass parameters to ecf model
Lambda = str.Lambda;
alpha = str.alpha;
th = str.th;
Th = str.Th;
Lh = str.Lh;
[betap, ferror] = m2mor(DA);
[np, mp, pm1] = size(Lambda);
p = pm1 - 1;
[na, ma] = size(alpha);
[nt, mt, qm1] = size(th);
q = qm1 - 1;
[nT, mT, Qm1] = size(Th);
Q = Qm1 - 1;
nlh = length(Lh);
np2 = np * np;
j = 0;
for i = 1:p
    Lambda(:, :, i+1) = reshape(x((j + i - 1)*np2+1:(j + i)*np2), np, np);
end
j = j + p * np2;
nama = na * ma;
alpha = reshape(x(j+1:j+nama), na, ma);
j = j + nama;
for i = 1:q
    th(:, :, i+1) = reshape(x(j+(i - 1)*np2+1:j+i*np2), np, np);
end
j = j + q * np2;
for i = 1:Q
    Th(:, :, i+1) = reshape(x(j+(i - 1)*np2+1:j+i*np2), np, np);
end
% phi,th,Phi,Th
Lh(2:end) = x(end-nlh+2:end);
L = zeros(np);
j = 0;
for i = 1:np
    inp = i:np;
    m = length(inp);
    L(inp, i) = Lh(j+1:j+m);
    j = j + m;
end
% L
