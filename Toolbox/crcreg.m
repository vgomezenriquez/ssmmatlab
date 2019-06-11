function [nr1, ns1, nr, ns] = crcreg(y, s, maxr)
%
% This function applies the CRC criterion to the y series. It is based on
% the paper "A Strongly Consistent Criterion to Decide Between I(1) and
% I(0) Processes Based on Different Convergence Rates" by V?ctor G?mez,
% (2013), Communications in Statistics - Simulation and Computation, 42,
% pp. 1848-1864.
%
% Input arguments:
% y    : series
% s    : number of seasons
% maxr : maximum regular differencing order considered
%
%
% Output arguments:
% nr1: number of regular differences found in the first step
% ns1: number of seasonal differences found in the first step
% nr: number of regular differences found
% ns: number of seasonal differences found
%
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

NN = length(y);
ols = 0;
a = 2;
if s <= 1, s = 0;
end
n = NN;
can = .13; %.11
if (s == 0)
    pord = 6;
    alphamax = .499;
    alpha1 = min(alphamax, .5-1/n);
    alpha2 = min(alphamax, .5-1/(n^.55));
else
    pord = max(6, floor(s/2)+1);
    alphamax = .499;
    alpha1 = min(alphamax, .5-1/(n^.95)); %.55 .65 .95
    alpha2 = min(alphamax, .5-1/(n^.90)); %.28  .51 .75 .90
    alpha2r = min(alphamax, .5-1/(n^.38)); %.33  .38 .47
end
%the case s=2 and maxr=2 would require special treatment in the first stage
maxs = 0;
if s > 1
    maxs = 1;
    p = min(pord, s-1);
else
    p = pord;
end
nr = 0;
ns = 0;
Y = ones(NN, 1);
%
% first step: AR model
%
q = 0;
ps = maxs;
qs = 0;
flag = 1;
hm = max(n^(-alphamax), n^(-alpha1));
hm1 = 1 - hm;
while flag == 1
    flag = 0;
    %  beta=bols(y,Y);
    beta = mean(y);
    x = inest([y, Y], beta, s, 0, p, ps, q, qs, 0, ols, a);
    rts = roots([1, x(1:p)]);
    if (ns < maxs)
        if (-x(p+1) > hm1)
            ns = ns + 1;
            NN = NN - s;
            y = diferm(y, s);
            Y = ones(NN, 1);
            flag = 1;
        end
    end
    if (flag == 0) && (nr < maxr)
        for j = 1:p
            if (real(rts(j)) > hm1) && (abs(imag(rts(j))) < hm)
                nr = nr + 1;
                NN = NN - 1;
                y = diferm(y, 1);
                Y = ones(NN, 1);
                flag = 1;
                break;
            end
        end
    end
end
%
% degrees of differencing found in the first step
%
nr1 = nr;
ns1 = ns;
if nr1 == maxr && ns1 == maxs, return;
end
%
% second step: ARMA model
%
p = 1;
q = 1;
ps = 0;
qs = 0;
flag = 0;
if ((nr < maxr) || (ns < maxs))
    flag = 1;
end
if (s > 1)
    ps = 1;
    qs = 1;
else
    s = 0;
end
pps = p + ps;
ppsq = pps + q;
hms = max(n^(-alphamax), n^(-alpha2));
hm1s = 1 - hms;
if (s > 1)
    hmsr = max(n^(-alphamax), n^(-alpha2r));
    hm1sr = 1 - hmsr;
end
while flag == 1
    flag = 0;
    yc = y - Y * mean(y);
    [x, laphi, x2] = hanris(yc, s, 0, p, ps, q, qs, 0, ols, a);
    chkar = chkroots(x(1:p+ps), p, ps, 0, 0, 0); %check stationarity
    if chkar == 1
        x = x2; %second step estimator in HR method, otherwise third step estimator
    end
    Rr = -x(p);
    %
    % increase differencing one by one, starting with seasonal differencing
    %
    if s > 1
        Rs = -x(p+ps);
        %   ns,nr
        %   x,Rs,hm1s,hm1sr,abs(x(pps)-x(ppsq+qs))
        if (ns < maxs) && ((Rs > hm1s) && (abs(x(pps)-x(ppsq+qs)) > can))
            ns = ns + 1;
            NN = NN - s;
            y = diferm(y, s);
            Y = ones(NN, 1);
            if (nr < maxr)
                flag = 1;
            end
        end
        if (nr < maxr) && (flag == 0) && ((Rr > hm1sr) && (abs(x(p)-x(ppsq)) > can))
            nr = nr + 1;
            NN = NN - 1;
            y = diferm(y, 1);
            Y = ones(NN, 1);
            if ((nr < maxr) || (ns < maxs))
                flag = 1;
            end
        end
    elseif (nr < maxr) && (((Rr > hm1s) && (abs(x(p)-x(ppsq))) > can))
        nr = nr + 1;
        NN = NN - 1;
        y = diferm(y, 1);
        Y = ones(NN, 1);
        if (nr < maxr)
            flag = 1;
        end
    end
end
