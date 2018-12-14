function Y = genfixseaspat(modescr, n)
%************************************************************************
% This function creates regression variables corresponding to fixed
% seasonal patterns of the form
%
%  s_t = a*cos(w*t) + b*sin(w*t),
%
% where w=2*pi*k/n
%
%    INPUTS:
%  modescr : structure with the following fields:
%             .seas : number of seasonal patterns
%            .seasp : cell array containing the pairs [per_j,m_j]
%                     for the seasonal patterns, where per_j is the period
%                     and m_j is the number of harmonics in the j-th
%                     seasonal pattern
%        n : desired length for the regression variables
%-------------------------------------------------------------------------
%
%   OUTPUTS:
%   Y : matrix with regression variables
%
%*************************************************************************
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
%*************************************************************************

if ~isfield(modescr, 'seas')
    seas = 0;
else
    seas = modescr.seas;
end
if ~isfield(modescr, 'seasp')
    seasp = [];
else
    seasp = modescr.seasp;
end


ncol = 0;
if seas > 0 %trigonometric seasonality
    %loop over the number of seasonal patterns
    for i = 1:seas
        per = seasp{i}(1);
        nh = seasp{i}(2);
        for j = 1:nh - 1
            ncol = ncol + 2;
        end
        if (mod(per, 2) == 0) && (nh == floor(per/2))
            ncol = ncol + 1;
        else
            ncol = ncol + 2;
        end
    end
else
    Y = [];
    return
end

if (ncol > 0)
    Y = zeros(n, ncol);
else
    Y = [];
    return
end


ncol = 0;
%loop over the number of seasonal patterns
for i = 1:seas
    per = seasp{i}(1);
    nh = seasp{i}(2);
    for j = 1:nh - 1
        omega = 2 * pi * j / per;
        for k = 1:n
            Y(k, ncol+1) = cos(omega*k);
            Y(k, ncol+2) = sin(omega*k);
        end
        ncol = ncol + 2;
    end
    if (mod(per, 2) == 0) && (nh == floor(per/2))
        for k = 1:n
            Y(k, ncol+1) = cos(pi*k);
        end
    else
        omega = 2 * pi * nh / per;
        for k = 1:n
            Y(k, ncol+1) = cos(omega*k);
            Y(k, ncol+2) = sin(omega*k);
        end
        ncol = ncol + 2;
    end
end
