function z = rootsarma(x, parm)
%
% given an ARMA model, this function computes a matrix containing the roots
% of the AR and MA polynomials, their arguments and their periods.
%
% Input arguments:
% x    : array containing model parameters
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
% .pvar:  array containing the indices of variable parameters
% .pfix:  array containing the indices of fixed parameters
% .ninput: number of inputs
% .delay: array with the delays of the input filters
% .ma: array with the ma parameters of the input filters
% .ar: array with the ar parameters of the input filters
%
% Output arguments:
% z: an array with four columns containing the real and imaginary parts of
%    the roots, their arguments and their periods
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

z = [];
tpi = 2 * pi;
pps = p + ps;
ppsq = pps + q;
ppsqqs = ppsq + qs;
ppsqqsqS = ppsqqs + qS;
if p > 0
    rts = roots([1, x(1:p)]);
    rrts = real(rts);
    irts = imag(rts);
    arts = abs(rts);
    argu = atan(irts./rrts);
    per = Inf(p, 1);
    for i = 1:p
        if argu(i) ~= 0.0
            per(i) = tpi / argu(i);
        end
    end
    z = [z; rrts, irts, arts, argu, per];
end
if ps > 0
    rts = roots([1, x(p+1:pps)]);
    rrts = real(rts);
    irts = imag(rts);
    arts = abs(rts);
    argu = atan(irts./rrts);
    per = Inf(ps, 1);
    for i = 1:ps
        if argu(i) ~= 0.0
            per(i) = tpi / argu(i);
        end
    end
    z = [z; rrts, irts, arts, argu, per];
end
if q > 0
    rts = roots([1, x(pps+1:ppsq)]);
    rrts = real(rts);
    irts = imag(rts);
    arts = abs(rts);
    argu = atan(irts./rrts);
    per = Inf(q, 1);
    for i = 1:q
        if argu(i) ~= 0.0
            per(i) = tpi / argu(i);
        end
    end
    z = [z; rrts, irts, arts, argu, per];
end
if qs > 0
    rts = roots([1, x(ppsq+1:ppsqqs)]);
    rrts = real(rts);
    irts = imag(rts);
    arts = abs(rts);
    argu = atan(irts./rrts);
    per = Inf(qs, 1);
    for i = 1:qs
        if argu(i) ~= 0.0
            per(i) = tpi / argu(i);
        end
    end
    z = [z; rrts, irts, arts, argu, per];
end
if qS > 0
    rts = roots([1, x(ppsqqs+1:ppsqqsqS)]);
    rrts = real(rts);
    irts = imag(rts);
    arts = abs(rts);
    argu = atan(irts./rrts);
    per = Inf(qS, 1);
    for i = 1:qS
        if argu(i) ~= 0.0
            per(i) = tpi / argu(i);
        end
    end
    z = [z; rrts, irts, arts, argu, per];
end
if ninput > 0
    cont = ppsqqsqS;
    for i = 1:ninput
        rts = roots([x(cont+1:cont+1+parm.ma(i))]);
        cont = cont + 1 + parm.ma(i);
        rrts = real(rts);
        irts = imag(rts);
        arts = abs(rts);
        nel = length(rrts);
        argu = atan(irts./rrts);
        per = Inf(nel, 1);
        for j = 1:parm.ma(i)
            if argu(j) ~= 0.0
                per(j) = tpi / argu(j);
            end
        end
        z = [z; rrts, irts, arts, argu, per];
        rts = roots([1, x(cont+1:cont+parm.ar(i))]);
        cont = cont + parm.ar(i);
        rrts = real(rts);
        irts = imag(rts);
        arts = abs(rts);
        nel = length(rrts);
        argu = atan(irts./rrts);
        per = Inf(nel, 1);
        for j = 1:parm.ar(i)
            if argu(j) ~= 0.0
                per(j) = tpi / argu(j);
            end
        end
        z = [z; rrts, irts, arts, argu, per];
    end
end
