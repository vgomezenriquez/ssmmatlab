function iout = miout(omet, C, delta, schr, nrout, nind, tip, Yo, ornames)
%
% function to create a structure containing parameters about outlier
% detection
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

iout.omet = omet;
iout.C = C;
iout.delta = delta;
iout.schr = schr;
iout.nrout = nrout;
iout.nind = nind;
iout.tip = tip;
iout.Yo = Yo;
iout.ornames = ornames;