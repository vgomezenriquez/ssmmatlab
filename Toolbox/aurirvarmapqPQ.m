function [str, ferror] = aurirvarmapqPQ(str, nr, DA)
% PURPOSE: given a structure containing information about a VARMA model, it
% adds the unit roots information given by nr and DA
%---------------------------------------------------
% USAGE: str = fixvarmaxpqPQ(str)
% where:str  = a structure created with function suvarmapqPQ
%        nr  = number of unit roots in the model
%        DA  = matrix of the form [DAr Indxr], where DAr is the
%              parameterization of betaor (the unit root part), and Indxr
%              is an index vector to identify the l.i. rows of DAr.
%---------------------------------------------------
% RETURNS: str = a structure containing model information
%          where the following fields have been added:
%         .nr  : number of unit roots in the model
%         .ns  : number of seasonal unit roots in the model (not used)
%         .xd  : parameter vector for betaor, including fixed and variable
%                parameters
%       .nparmd: number of parameters for the unit root part
%         .xid : array of ones and zeros, as in xi, for the unit root part
%         .xvd : array of parameters to estimate for the unit root part
%         .xfd : array of fixed parameters for the unit root part
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

str.nr = nr;
str.DA = DA;
str.ns = 0;

%obtain parameters in differencing polynomial
[nd, md] = size(DA);
xd = [];
DAn = DA;
my = nd;
mds = md;
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
