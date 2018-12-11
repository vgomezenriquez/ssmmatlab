function [str, ferror] = aurivarmapqPQ(str, nr, ns, DA)
% PURPOSE: given a structure containing information about a VARMA model, it
% adds the unit roots information given by nr, ns and DA ( the output from
% mcrcreg)
%---------------------------------------------------
% USAGE: str = fixvarmaxpqPQ(str)
% where:   str      = a structure created with function suvarmapqPQ
%---------------------------------------------------
% RETURNS: str = a structure containing model information
%---------------------------------------------------
% Copyright (c) 21 July 2003 by Victor Gomez
% Ministerio de Hacienda, Direccion Gral. de Presupuestos,
% Subdireccion Gral. de Analisis y P.E.,
% Alberto Alcocer 2, 1-P, D-34, 28046, Madrid, SPAIN.
% Phone : +34-915835439
% E-mail: vgomez@sgpg.meh.es
%
% The author assumes no responsibility for errors or damage resulting from the use
% of this code. Usage of this code in applications and/or alterations of it should
% be referenced. This code may be redistributed if nothing has been added or
% removed and no money is charged. Positive or negative feedback would be appreciated.
%

ferror = 0;
% if (nr+ns == 0)
%  return
% end

str.nr = nr;
str.ns = ns;
str.DA = DA;

%obtain parameters in differencing polynomial
[nd, md] = size(DA);
xd = [];
DAn = DA;
my = nd;
mds = md;
if (ns > 0)
    Indxs = DA(:, md);
    sindxs = sum(Indxs);
    sncols = my - sindxs;
    mdm1 = md - 1;
    mds = mdm1 - sncols;
    mind = mdm1 - sncols + 1:mdm1;
    xd = zeros(1, sncols*sindxs);
    cont = 0;
    for i = 1:my
        if (Indxs(i) == 1)
            xd(cont*sncols+1:(cont + 1)*sncols) = DA(i, mind);
            DAn(i, mind) = NaN(size(DA(i, mind)));
            cont = cont + 1;
        end
    end
end
if (nr > 0)
    Indxr = DA(:, mds);
    sindxr = sum(Indxr);
    rncols = my - sindxr;
    mdsm1 = mds - 1;
    mind = mdsm1 - rncols + 1:mdsm1;
    lxs = length(xd);
    xd = [xd, zeros(1, rncols*sindxr)];
    cont = 0;
    for i = 1:my
        if (Indxr(i) == 1)
            xd(lxs+cont*rncols+1:lxs+(cont + 1)*rncols) = DA(i, mind);
            DAn(i, mind) = NaN(size(DA(i, mind)));
            cont = cont + 1;
        end
    end
end
% xd
nparmd = length(xd);
str.xd = xd;
str.nparmd = nparmd;
str.DAn = DAn;
xid = ones(size(xd));
xvd = xd;
xfd = [];
str.xid = xid;
str.xvd = xvd;
str.xfd = xfd;
