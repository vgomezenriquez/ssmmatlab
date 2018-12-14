function oparm = leapid(y, Y, infm, parm, ser, ols, a, fid, fmarqdt)
%
% this function identifies the leap year period using BIC
%
% Copyright (c) 21 July 2003 by Victor Gomez
% Ministerio de Hacienda, Direccion Gral. de Presupuestos,
% Subdireccion Gral. de Analisis y P.E.,
% Alberto Alcocer 2, 1-P, D-34, 28046, Madrid, SPAIN.
% Phone : +34-915835439
% E-mail: vgomez@sepg.minhac.es
%
% The author assumes no responsibility for errors or damage resulting from the use
% of this code. Usage of this code in applications and/or alterations of it should
% be referenced. This code may be redistributed if nothing has been added or
% removed and no money is charged. Positive or negative feedback would be appreciated.
%

if ~isstruct(parm)
    error('leapid: requires a parameter structure');
end;
oparm = parm;
ny = length(y);
s = parm.s;
S = parm.S;
dr = parm.dr;
ds = parm.ds;
dS = parm.dS;
pvar = parm.pvar;
pfix = parm.pfix;
p = oparm.p;
q = oparm.q;
ps = oparm.ps;
qs = oparm.qs;
qS = oparm.qS;
bg_year = ser.bg_year;
bg_per = ser.bg_per;
freq = ser.freq;
bicm = 1.d10;
leapval = [0, 1];
for i = 1:2
    leap = leapval(i);
    if leap > 0, Yle = genleap(bg_year, bg_per, ny, freq);
        YY = [Y, Yle];
    else, YY = Y;
    end
    est = 1;
    x0 = cinest(y, YY, parm, est, ols, a, 0, fid);
    xv = x0(pvar);
    xf = x0(pfix);
    if ~isempty(xv)
        xv = arimaopt(fmarqdt, fid, x0, xv, xf, y, YY, parm, infm, 0);
        x0(pvar) = xv;
        %  [F,e,beta,M]=residual2(x0,y,YY,s,dr,ds,p,ps,q,qs,1); nyd=ny-dr-ds*s;
    end
    [F, e, beta, M] = residual2x(x0, y, YY, s, S, dr, ds, dS, p, ps, q, qs, qS);
    nyd = ny - dr - ds * s - dS * S;
    nbeta = length(beta);
    dn = double(nyd);
    var = e' * e / dn;
    ldn = log(dn);
    bic = log(var) + double(nbeta) * ldn / dn;
    if (bic < bicm)
        bicm = bic;
        oparm.leap = leap;
    end
end
