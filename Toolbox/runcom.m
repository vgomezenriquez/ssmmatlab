function [n0, n1, nr, Tval] = runcom(X, N, Xmed)
%*************************************************************************
%This function computes the elements of a test of randomness based on runs.
%
%   INPUTS:
%       X : data vector
%       N : length of X
%    Xmed : median of X
%
%  OUTPUTS:
%      n0 : number of values of X smaller than Xmed (-)
%      n1 : number of values of X grater or equal to Xmed (+)
%      nr : total number of runs, (-) and (+).
%    Tval : approximate t-value of nr
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
%*************************************************************************

nr = 1;
l = 0;
if (X(1) >= Xmed), l = 1;
end
n1 = l;
n0 = 1 - l;
for i = 2:N
    ll = 0;
    if (X(i) >= Xmed), ll = 1;
    end
    n1 = n1 + ll;
    n0 = n0 + 1 - ll;
    if (ll ~= l)
        l = ll;
        nr = nr + 1;
    end
end
runm = 1.d0 + double(2*n1*n0) / double(n1+n0);
runstd = double(2*n1*n0) * double(2*n1*n0-n1-n0);
runstd = runstd / (double((n1 + n0)^2) * double(n1+n0-1));
runstd = sqrt(runstd);
if (runstd < 1.d-9), runstd = 1.d-9;
end
Tval = (double(nr) - runm) / runstd;
