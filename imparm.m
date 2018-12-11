function [s, S, dr, ds, dS, p, ps, q, qs, qS, ny, nreg, pfix, pvar] = imparm(parm)
%
% this function obtains the paramters in structure parm
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

s = parm.s;
S = parm.S;
dr = parm.dr;
ds = parm.ds;
dS = parm.dS;
p = parm.p;
ps = parm.ps;
q = parm.q;
qs = parm.qs;
qS = parm.qS;
ny = parm.ny;
nreg = parm.nreg;
pfix = parm.pfix;
pvar = parm.pvar;
