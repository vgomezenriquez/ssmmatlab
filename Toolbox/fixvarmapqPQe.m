function [str, ferror] = fixvarmapqPQe(str)
% PURPOSE: given a structure containing information about a VARMA model, it
% fixes the parameters in phi, Phi, th, Th and Lh that correspond to those
% elements in phin, Phin, thn, Thn and Lhn that are not equal to NaN.
%---------------------------------------------------
% USAGE: str = fixvarmaxpqPQe(str)
% where:   str      = a structure created with function suvarmapqPQ
%---------------------------------------------------
% RETURNS: str = a structure containing model information
%---------------------------------------------------
% Copyright (c) 21 July 2003 by Victor Gomez
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

ferror = 0;

if (~isfield(str, 'Lambdan') || ~isfield(str, 'alphan') || ...
        ~isfield(str, 'thn') || ~isfield(str, 'Thn') || ...
        ~isfield(str, 'Lhn') || ~isfield(str, 'DAn'))
    ferror = 1;
    disp('some of the fields Lambdan, alphan, thn, Thn, Lhn')
    disp('or DAn are not in fixvarmapqPQe');
    return
end

Lambdan = str.Lambdan;
thn = str.thn;
alphan = str.alphan;
Thn = str.Thn;
Lhn = str.Lhn;

[np, mp, pm1] = size(Lambdan);
[nt, mt, qm1] = size(thn);
[nT, mT, Qm1] = size(Thn);

%create vector for fixed model parameters
x = [];
for i = 2:pm1
    x = [x, vec(Lambdan(:, :, i))'];
end
x = [x, vec(alphan)'];
for i = 2:qm1
    x = [x, vec(thn(:, :, i))'];
end
for i = 2:Qm1
    x = [x, vec(Thn(:, :, i))'];
end
x = [x, Lhn(2:end)'];

nparm = str.nparm;
xi = str.xi;
xx = str.x;
xv = [];
xf = [];
for i = 1:nparm
    if isnan(x(i))
        xv = [xv, xx(i)];
    else
        xf = [xf, xx(i)];
        xi(i) = 0;
    end
end
str.xi = xi;
str.xv = xv;
str.xf = xf;
%add unit root information
DAn = str.DAn;
nr = str.nr;
%create vector for fixed model parameters
[nd, md] = size(DAn);
xd = [];
my = nd;
mds = md;
if (nr > 0)
    Indxr = DAn(:, mds);
    sindxr = sum(Indxr);
    rncols = my - sindxr;
    mdsm1 = mds - 1;
    mind = mdsm1 - rncols + 1:mdsm1;
    lxs = length(xd);
    xd = [xd, zeros(1, rncols*sindxr)];
    cont = 0;
    for i = 1:my
        if (Indxr(i) == 1)
            xd(lxs+cont*rncols+1:lxs+(cont + 1)*rncols) = DAn(i, mind);
            %    DAn(i,mind)=NaN(size(DAn(i,mind)));
            cont = cont + 1;
        end
    end
end
% xd

nparmd = length(xd);
xid = str.xid;
xxd = str.xd;
xvd = [];
xfd = [];
for i = 1:nparmd
    if isnan(xd(i))
        xvd = [xvd, xxd(i)];
    else
        xfd = [xfd, xxd(i)];
        xid(i) = 0;
    end
end
str.xid = xid;
str.xvd = xvd;
str.xfd = xfd;
