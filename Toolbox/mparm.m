function parm = mparm(s, S, dr, ds, dS, p, ps, q, qs, qS, ny, nreg, pfix, pvar, lam, flagm, ...
    trad, leap, east, dur, ninput, nmiss)
%
% function to create a structure containing ARIMA parameters
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

parm.s = s;
parm.S = S;
parm.dr = dr;
parm.ds = ds;
parm.dS = dS;
parm.p = p;
parm.ps = ps;
parm.q = q;
parm.qs = qs;
parm.qS = qS;
parm.ny = ny;
parm.nreg = nreg;
parm.pfix = pfix;
parm.pvar = pvar;
parm.lam = lam;
parm.flagm = flagm;
parm.trad = trad;
parm.leap = leap;
parm.east = east;
parm.dur = dur;
parm.ninput = ninput;
parm.inputv = [];
parm.delay = [];
parm.ma = [];
parm.ar = [];
parm.nmiss = nmiss;