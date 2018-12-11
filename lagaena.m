function [nlag, aenames] = lagaena(parm)
%
% function to create names for ARIMA estimation printing
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
s = parm.s;
ninput = parm.ninput;
S = parm.S;

nlag = [];
aenames = strvcat('Parameter');
if p > 0
    nlag = [nlag; ones(p, 1)];
    for i = 1:p, aenames = strvcat(aenames, ['phi', num2str(i)]);
    end
end
if ps > 0
    nlag = [nlag; s * ones(ps, 1)];
    for i = 1:ps, aenames = strvcat(aenames, ['phi', num2str(s*i)]);
    end
end
if q > 0
    nlag = [nlag; ones(q, 1)];
    for i = 1:q, aenames = strvcat(aenames, ['ma', num2str(i)]);
    end
end
if qs > 0
    nlag = [nlag; s * ones(qs, 1)];
    for i = 1:qs, aenames = strvcat(aenames, ['ma', num2str(s*i)]);
    end
end
if qS > 0
    nlag = [nlag; S * ones(qS, 1)];
    for i = 1:qS, aenames = strvcat(aenames, ['ma', num2str(S*i)]);
    end
end
if ninput > 0
    idx = parm.inputidx;
    for i = 1:ninput
        for j = parm.delay(i):parm.delay(i) + parm.ma(i)
            aenames = strvcat(aenames, ['omg', num2str(idx(i)), num2str(j)]);
        end
        for j = 1:parm.ar(i)
            aenames = strvcat(aenames, ['del', num2str(idx(i)), num2str(j)]);
        end
        eval(['nlag=[nlag; ones(parm.ma(', num2str(i), ')+1,1)];'])
        eval(['nlag=[nlag; ones(parm.ar(', num2str(i), '),1)];'])
    end
end
