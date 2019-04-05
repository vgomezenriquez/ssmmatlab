function oparm = armaid(y, parm, ols, a, maxpq, maxPQ)
%
% this function automatically identifies an ARMA model for a stationary
% series using the BIC criterion.
%
% Input arguments:
% y: vector containing the data
% parm: astructure containing model information, where
%   s:  seasonality
%   S:  second seasonality
% .dr: order of regular differencing
% .ds: order of differencing of order s
% .dS: order of differencing of order S
% .p:  initial AR order
% .ps: initial order of the AR of order s
% .q:  initial order of the regular MA
% .qs: initial order of the MA of order s (1 at most)
% .qS: initial order of the MA of order S (1 at most)
% ols: = 0, use the Levinson-Durbin algorithm in the Hannan-Rissanen method
%      = 1, use OLS in  the Hannan-Rissanen method
%   a: the exponent in log(n)^a for the length of the long AR in the
%      Hannan-Rissanen method. By default, a = 1.5.
% maxpq: the maximum orders of the regualr AR and MA polynomials
% maxpq: the maximum orders of the seasonal AR and MA polynomials
%
% Output arguments:
% oparm: a structure containing the same fields as par plus
% .p:  AR order
% .ps: order of the AR of order s
% .q:  order of the regular MA
% .qs: order of the MA of order s (1 at most)
% .pvar:  array containing the indices of variable parameters
% .pfix:  array containing the indices of fixed parameters
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

if ~isstruct(parm)
    error('armaid: requires a parameter structure');
end
C = .005;
oparm = parm;
s = parm.s;
S = parm.S;
dr = parm.dr;
ds = parm.ds;
dS = parm.dS;
est = 0;
Y = [];
qS = parm.qS;
[yd, beta] = diffest(y, Y, s, S, dr, ds, dS, est);
yd = yd - mean(yd) * ones(size(yd));
bicm = 1.d10;
sbic = bicm;
sparm.p = 0;
sparm.ps = 0;
sparm.q = 0;
sparm.qs = 0;
if s > 1
    p = maxpq;
    q = 0;
    sparm.p = oparm.p;
    sparm.q = oparm.q;
    oparm.p = p;
    oparm.q = q;
    for ps = 0:maxPQ
        for qs = 0:maxPQ
            [bic, ~] = updbic(yd, beta, s, S, p, ps, q, qs, qS, ols, a, bicm, oparm);
            if bic < bicm
                sbic = bicm;
                bicm = bic;
                sparm.ps = oparm.ps;
                sparm.qs = oparm.qs;
                oparm.ps = ps;
                oparm.qs = qs;
            end
        end
    end
    ps = oparm.ps;
    qs = oparm.qs;
    bicm = 1.d10;
    for p = 0:maxpq
        for q = 0:maxpq
            [bic, ~] = updbic(yd, beta, s, S, p, ps, q, qs, qS, ols, a, bicm, oparm);
            if bic < bicm
                sbic = bicm;
                bicm = bic;
                sparm.p = oparm.p;
                sparm.q = oparm.q;
                oparm.p = p;
                oparm.q = q;
                %     p,q,bic,pause
            end
        end
    end
    p = oparm.p;
    q = oparm.q;
    sp = sparm.p;
    sq = sparm.q;
    rbic = sbic - bicm;
    % p,ps,q,qs,sparm
%     bicm,sbic,rbic
    if rbic < C
        if ((sp + sq < p + q) && (sp + sq > 0)) 
            oparm.q = sq;
            oparm.p = sp;
        end
    end
    p = oparm.p;
    q = oparm.q;
    bicm = 1.d10;
    for ps = 0:maxPQ
        for qs = 0:maxPQ
            [bic, ~] = updbic(yd, beta, s, S, p, ps, q, qs, qS, ols, a, bicm, oparm);
            if bic < bicm
                sbic = bicm;
                bicm = bic;
                sparm.ps = oparm.ps;
                sparm.qs = oparm.qs;
                oparm.ps = ps;
                oparm.qs = qs;
            end
        end
    end
    ps = oparm.ps;
    qs = oparm.qs;
    sps = sparm.ps;
    sqs = sparm.qs;
    rbic = sbic - bicm;
    % p,ps,q,qs,sparm
%     bicm,sbic,rbic
    if rbic < C
        if  ((sps + sqs < ps + qs) && (sps + sqs > 0))
            oparm.qs = sqs;
            oparm.ps = sps;
        end
    end
else
    ps = 0;
    qs = 0;
    for p = 0:maxpq
        for q = 0:maxpq
            [bic, ~] = updbic(yd, beta, s, S, p, ps, q, qs, qS, ols, a, bicm, oparm);
            if bic < bicm
                sbic = bicm;
                bicm = bic;
                sparm.p = oparm.p;
                sparm.q = oparm.q;
                oparm.p = p;
                oparm.q = q;
            end
        end
    end
    p = oparm.p;
    q = oparm.q;
    sp = sparm.p;
    sq = sparm.q;
    rbic = sbic - bicm;
    % p,ps,q,qs,sparm
    % bicm,sbic,rbic
    if rbic < C
        if ((sp + sq < p + q) && (sp + sq > 0)) 
            oparm.q = sq;
            oparm.p = sp;
        end
    end
end
p = oparm.p;
ps = oparm.ps;
q = oparm.q;
qs = oparm.qs;
if p + q + ps + qs == 0 %added on 25-5-2006
    oparm.q = 1;
    q = 1;
end
nr = p + ps + q + qs + qS; %number of arma parameters
pfix = [];
pvar = [1:nr]; %fixed and free parameters
oparm.pfix = pfix;
oparm.pvar = pvar;
