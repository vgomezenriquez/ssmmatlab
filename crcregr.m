function [nr1, nr] = crcregr(y, maxr)
%
% This function applies the CRC criterion to the y series. It is based on
% the paper "A Strongly Consistent Criterion to Decide Between I(1) and
% I(0) Processes Based on Different Convergence Rates" by Víctor Gómez,
% (2013), Communications in Statistics - Simulation and Computation, 42,
% pp. 1848-1864.
%
% Input arguments:
% y: series
% maxr: maximum regular differencing order considered
%
% Output arguments:
% nr1: number of regular differences found in the first step
% nr: number of regular differences found
%
%
% Copyright (c) 21 July 2014 by Victor Gomez
% Ministerio de Hacienda y A.P., Direccion Gral. de Presupuestos,
% Subdireccion Gral. de Analisis y P.E.,
% Alberto Alcocer 2, 1-P, D-34, 28046, Madrid, SPAIN.
% Phone : +34-915835439
% E-mail: vgomez@sepg.minhap.es
%
% The author assumes no responsibility for errors or damage resulting from the use
% of this code. Usage of this code in applications and/or alterations of it should
% be referenced. This code may be redistributed if nothing has been added or
% removed and no money is charged. Positive or negative feedback would be appreciated.
%


NN = length(y);
ols = 0;
a = 2;
n = NN;
can = .13;
pord = 6; %.11
alphamax = .499;
alpha1 = min(alphamax, .5-1/n);
alpha2 = min(alphamax, .5-1/(n^.55));
p = pord;
nr = 0;
Y = ones(NN, 1);
%
% first step: AR model
%
q = 0;
ps = 0;
qs = 0;
flag = 1;
s = 0;
hm = max(n^(-alphamax), n^(-alpha1));
hm1 = 1 - hm;
while flag == 1
    flag = 0;
    %  beta=bols(y,Y);
    beta = mean(y);
    x = inest([y, Y], beta, s, 0, p, ps, q, qs, 0, ols, a);
    rts = roots([1, x(1:p)]);
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
if nr1 == maxr
    return;
end
%
% second step: ARMA model
%
p = 1;
q = 1;
ps = 0;
qs = 0;
flag = 0;
if (nr < maxr)
    flag = 1;
end
s = 0;
pps = p + ps;
ppsq = pps + q;
hms = max(n^(-alphamax), n^(-alpha2));
hm1s = 1 - hms;
while flag == 1
    flag = 0;
    yc = y - Y * mean(y);
    [x, laphi, x2] = hanris(yc, s, 0, p, ps, q, qs, 0, ols, a);
    chkar = chkroots(x(1:p+ps), p, ps, 0, 0, 0); %check stationarity
    if chkar == 1
        x = x2;
    end
    Rr = -x(p);
    if (nr < maxr) && ((Rr > hm1s) && (abs(x(p)-x(ppsq))) > can)
        nr = nr + 1;
        NN = NN - 1;
        y = diferm(y, 1);
        Y = ones(NN, 1);
        if (nr < maxr), flag = 1;
        end
    end
end
