function [yd, xvv, xff, DA, Dr, Ds, ferror] = pr2varmapqPQd(y, xv, xf, str)
% PURPOSE: given a structure containing information about a VARMA model
% with unit roots, it passes the parameters to the state space form.
%---------------------------------------------------
% USAGE: [yd,xvv,xff,DA,Dr,Ds,ferror] = pr2varmapqPQd(y,xv,xf,str)
% where:   y        = an (n x neqs) matrix containing the data
%          xv       = a vector containing the parameters to be estimated
%          xf       = a vector containing the fixed parameters
%                     be estimated, =0, not
%          str      = a structure containing the initial model information
%---------------------------------------------------
%---------------------------------------------------
% where:   yd     = the differenced series
%          xvv    = the variable paramaters for VARMA model
%          xff    = the fixed paramaters for VARMA model
%          DA     = the matrix containing the parameterization of the
%                   differencing polynomial
%          Dr     = the regular 'differencing' polynomial
%          Ds     = the seasonal 'differencing' polynomial (not used)
%          ferror   = flag for errors
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

[ny, my] = size(y);
fur = 0;
if isfield(str, 'nr')
    nr = str.nr;
    fur = 1;
end
if isfield(str, 'ns')
    ns = str.ns;
    if (fur == 0)
        disp('pr2varmapqPQd: if there are unit roots, fields nr, ns and DA should be present in input structure')
        ferror = 1;
        return
    end
    fur = 1;
elseif (fur == 1)
    disp('pr2varmapqPQd: if there are unit roots, fields nr, ns and DA should be present in input structure')
    ferror = 1;
    return
end
if isfield(str, 'DA')
    DA = str.DA;
    if (fur == 0)
        disp('pr2varmapqPQd: if there are unit roots, fields nr, ns and DA should be present in input structure')
        ferror = 1;
        return
    end
    fur = 1;
elseif (fur == 1)
    disp('pr2varmapqPQd: if there are unit roots, fields nr, ns and DA should be present in input structure')
    ferror = 1;
    return
end
if (fur == 1)
    if (~isfield(str, 'xd'))
        ferror = 1;
        disp('the field xd is not in pr2varmapqPQd');
        return
    end
    if (~isfield(str, 'xid'))
        ferror = 1;
        disp('the field xid is not in pr2varmapqPQd');
        return
    else
        xid = str.xid;
    end
    if (~isfield(str, 'nparmd'))
        ferror = 1;
        disp('the field nparmd is not in pr2varmapqPQd');
        return
    end
    if (~isfield(str, 'DAn'))
        ferror = 1;
        disp('the field DAn is not in pr2varmapqPQd');
        return
    end
    if (~isfield(str, 'xvd'))
        ferror = 1;
        disp('the field xvd is not in pr2varmapqPQd');
        return
    else
        xvd = str.xvd;
    end
    if (~isfield(str, 'xfd'))
        ferror = 1;
        disp('the field xfd is not in pr2varmapqPQd');
        return
    else
        xfd = str.xfd;
    end
    nxvd = length(xvd);
    nxv = length(xv);
    nxf = length(xf);
    nxfd = length(xfd);
    xvv = xv(nxvd+1:nxv);
    xff = xf(nxfd+1:nxf);
    xvd = xv(1:nxvd);
    xfd = xf(1:nxfd);
    if isempty(xvd)
        xd = xfd;
    else
        if isempty(xfd)
            xd = xvd;
        else
            xd = zeros(size(xid));
            xd(xid > 0) = xvd;
            xd(xid == 0) = xfd;
        end
    end
    [nd, md] = size(DA);
    mds = md;
    if (ns > 0)
        Indxs = DA(:, md);
        sindxs = sum(Indxs);
        sncols = my - sindxs;
        mdm1 = md - 1;
        mds = mdm1 - sncols;
        mind = mdm1 - sncols + 1:mdm1;
        cont = 0;
        for i = 1:my
            if (Indxs(i) == 1)
                DA(i, mind) = xd(cont*sncols+1:(cont + 1)*sncols);
                cont = cont + 1;
            end
        end
    else
        sncols = 0;
        sindxs = 0;
    end
    if (nr > 0)
        Indxr = DA(:, mds);
        sindxr = sum(Indxr);
        rncols = my - sindxr;
        mdsm1 = mds - 1;
        mind = mdsm1 - rncols + 1:mdsm1;
        if (ns > 0)
            lxs = sncols * sindxs;
        else
            lxs = 0;
        end
        cont = 0;
        for i = 1:my
            if (Indxr(i) == 1)
                DA(i, mind) = xd(lxs+cont*rncols+1:lxs+(cont + 1)*rncols);
                cont = cont + 1;
            end
        end
    end
    seas = str.freq;
    [yd, Dr, Ds, ferror] = param2mdp(y, DA, nr, ns, seas);
else
    xvv = xv;
    xff = xf;
    yd = y;
end
