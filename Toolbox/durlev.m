function [fi, pc] = durlev(c0, cv)
%*************************************************************************
% This function applies the Durbin-Levinson algorithm to fit an AR model of
% order p = length(cv), given the autocovariances cv and the variance c0.
% As a byproduct, it also gives the partial autocorrelation coefficients.
%
%  INPUTS:
%     c0 : variance
%     cv : an (1 x p) vector containing the i-th covariances, i=1,...,p.
%
% OUTPUTS:
%     fi : an (1 x p) vector containing the AR(p) polynomial
%     pc : an (1 x p) vector containing the partial correlation coefficients
%
% notation of fi is that of Box and Jenkins:
% y(t)-fi(1)*y(t-1)-...-fi(p)*y(t-p) = a(t)
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

p = length(cv);
fi = zeros(size(cv));
pc = fi;
fi(1) = cv(1) / c0;
pc(1) = fi(1);
v = (1 - fi(1)^2) * c0;
if p == 1
    return
else
    for i = 2:p
        fi(i) = (cv(i) - fi(1:i-1)' * flipud(cv(1:i-1))) / v;
        pc(i) = fi(i);
        v = (1 - fi(i)^2) * v;
        fi(1:i-1) = fi(1:i-1) - fi(i) * flipud(fi(1:i-1));
    end
end
