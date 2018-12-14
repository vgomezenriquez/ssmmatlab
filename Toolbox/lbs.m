function [qstat, pval, df, sea] = lbs(ne, p, r, nr)
%
% this function obtains the Q-statistics and their p-values for a sequence of integers p
%
% input arguments:
%              ne: length of the variable
%               p: a sequence of integers, i.e. p=1:lag
%               r: the covariance sequence
%              nr: number of parameters
% output arguments:
%           qstat: an array containing the Q-statistics
%            pval: an array containing the P-values
%              df: an array containing the degrees of freedom
%             sea: standard errors
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


if (nargin ~= 4) error('Wrong number of arguments to lbs');
end;
np = length(p);
qstat = zeros(np, 1);
pval = zeros(np, 1);
df = zeros(np, 1);
sea = zeros(np, 1);
sea(1) = 1 / sqrt(ne);
s = 0; %standard errors
for i = 1:np - 1
    s = s + r(i) * r(i);
    sea(i+1) = sqrt((1 + 2 * s)/ne);
end
for i = 1:np;
    qstat(i) = sum((r(1:i).^2)./[ne - 1:-1:ne - i]') * ne * (ne + 2);
    df(i) = max(0, i-nr); %degrees of freedom
    if df(i) > 0
        pval(i) = 1 - gammp(df(i)*.5, qstat(i)*.5);
    else
        pval(i) = 1;
    end
end;
