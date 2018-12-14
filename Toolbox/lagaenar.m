function aenames = lagaenar(parm)
%
% this function generates an string array containing the names for the
% roots of all the polynomials of a transfer function model. These include
% the ARMA polynomials and the rational input filters.
%
% Input arguments:
% parm : a structure where
% .s:  seasonality
% .S:  second seasonality
% .p:  AR order
% .ps: order of the AR of order s
% .q:  order of the regular MA
% .qs: order of the MA of order s (1 at most)
% .qS: order of the MA of order S (1 at most)
% .dr: order of regular differencing
% .ds: order of differencing of order s
% .dS: order of differencing of order S
% .ninput: number of inputs
% .inputidx: an index array containing the indices of the input variables
% .delay: array with the delays of the input filters
% .ma: array with the ma parameters of the input filters
% .ar: array with the ar parameters of the input filters
%
% Output arguments:
% aenames : a string array containing the names for the
% roots of all the polynomials in the transfer function model
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

p = parm.p;
q = parm.q;
ps = parm.ps;
qs = parm.qs;
qS = parm.qS;
ninput = parm.ninput;

aenames = strvcat('Root');
if p > 0
    for i = 1:p, aenames = strvcat(aenames, ['rarroot', num2str(i)]);
    end
end
if ps > 0
    for i = 1:ps, aenames = strvcat(aenames, ['sarroot', num2str(i)]);
    end
end
if q > 0
    for i = 1:q, aenames = strvcat(aenames, ['rmaroot', num2str(i)]);
    end
end
if qs > 0
    for i = 1:qs, aenames = strvcat(aenames, ['smaroot', num2str(i)]);
    end
end
if qS > 0
    for i = 1:qS, aenames = strvcat(aenames, ['Smaroot', num2str(i)]);
    end
end
if ninput > 0
    idx = parm.inputidx;
    for i = 1:ninput
        for j = 1:parm.ma(i)
            aenames = strvcat(aenames, ['omgroot', num2str(idx(i)), num2str(j)]);
        end
        for j = 1:parm.ar(i)
            aenames = strvcat(aenames, ['delroot', num2str(idx(i)), num2str(j)]);
        end
    end
end
