function str = restrcmodel(s, m, seas, ordersr, orderss)
% PURPOSE: given the regular and seasonal orders of a VARMAX(p,q,r)
%          (P,Q,R)_seas model, this function creates a structure containing
%          the polynomial and state space forms of the model with NaNs for
%          the parameters to be estimated and zero otherwise.
%---------------------------------------------------
% USAGE:  str = restrcmodel(s,m,seas,ordersr,orderss)
% where:    s       = an integer, the dimension of the output y_t
%           m       = an integer, the dimension of the input  x_t
%        seas       = seasonality
%        ordersr    = a 1 x 3 array containing the regular VARMAX orders
%        orderss    = a 1 x 3 array containing the seasonal VARMAX orders
%---------------------------------------------------
% RETURNS: str, a structure containing model information. The fields of
%          str are the same than those in the structure returned by
%          function matechelon, but with the zero restrictions of the
%          VARMAX model imposed. The VARMAX model is assumed to follow an
%          VARMAX model in echelon form with all Kronecker indices equal to
%          max(p,q,r) + seas*max(P,Q,R) and with the zero restrictions
%          imposed.
%---------------------------------------------------
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

%transform multiplicative to non-multiplicative model
maxordr = max(ordersr);
p = ordersr(1);
q = ordersr(2);
r = ordersr(3);
maxords = max(orderss);
if (seas > 1)
    %  maxord=maxordr + seas*maxords;
    P = orderss(1);
    Q = orderss(2);
    R = orderss(3);
    maxord = max([p + seas * P, q + seas * Q, r + seas * R]);
else
    maxord = maxordr;
    P = 0;
    Q = 0;
    R = 0;
end

kro = repmat(maxord, 1, s);
str = matechelon(kro, s, m);
s2 = s * s;
sm = s * m;

%restrictions in the model
%regular part
nparm = str.nparm;
if (maxordr + 1 >= p + 2)
    str.phi(:, :, p+2:maxordr+1) = zeros(s, s, maxordr-p);
    nparm = nparm - (maxordr - p) * s2;
end
if (maxordr + 1 >= q + 2)
    str.theta(:, :, q+2:maxordr+1) = zeros(s, s, maxordr-q);
    nparm = nparm - (maxordr - q) * s2;
end
if ~isempty(str.gamma)
    if (maxordr + 1 >= r + 2)
        str.gamma(:, :, r+2:maxordr+1) = zeros(s, m, maxordr-r);
        nparm = nparm - (maxordr - r) * sm;
    end
end
%seasonal part
if (seas > 1) && (maxords > 0)
    if (seas >= maxordr + 2)
        str.phi(:, :, maxordr+2:seas) = zeros(s, s, seas-maxordr-1);
        nparm = nparm - (seas - maxordr - 1) * s2;
    end
    if ((maxord + 1 >= seas + 1) && (P == 0))
        str.phi(:, :, seas+1:maxord+1) = zeros(s, s, maxord-seas+1);
        nparm = nparm - (maxord - seas + 1) * s2;
    end
    if ((maxord + 1 >= seas + p + 2) && (P > 0))
        str.phi(:, :, seas+p+2:maxord+1) = zeros(s, s, maxord-seas-p);
        nparm = nparm - (maxord - seas - p) * s2;
    end
    if (seas >= maxordr + 2)
        str.theta(:, :, maxordr+2:seas) = zeros(s, s, seas-maxordr-1);
        nparm = nparm - (seas - maxordr - 1) * s2;
    end
    if ((maxord + 1 >= seas + 1) && (Q == 0))
        str.theta(:, :, seas+1:maxord+1) = zeros(s, s, maxord-seas+1);
        nparm = nparm - (maxord - seas + 1) * s2;
    end
    if ((maxord + 1 >= seas + q + 2) && (Q > 0))
        str.theta(:, :, seas+q+2:maxord+1) = zeros(s, s, maxord-seas-q);
        nparm = nparm - (maxord - seas - q) * s2;
    end
    if ~isempty(str.gamma)
        if (seas >= maxordr + 2)
            str.gamma(:, :, maxordr+2:seas) = zeros(s, m, seas-maxordr-1);
            nparm = nparm - (seas - maxordr - 1) * sm;
        end
        if ((maxord + 1 >= seas + 1) && (R == 0))
            str.gamma(:, :, seas+1:maxord+1) = zeros(s, m, maxord-seas+1);
            nparm = nparm - (maxord - seas + 1) * sm;
        end
        if ((maxord + 1 >= seas + r + 2) && (R > 0))
            str.gamma(:, :, seas+r+2:maxord+1) = zeros(s, m, maxord-seas-r);
            nparm = nparm - (maxord - seas - r) * sm;
        end
    end
end
str.nparm = nparm;
