function [yd, Dr, Ds, ferror] = param2mdp(y, DA, nr, ns, seas)
%
%
% This function obtains the series
%         yd_t = D(B)*y_t,
% where D(z)=Dr(z)*Ds(z) is a polynomial matrix compatible
% with y_t and B is the backshift operator, By_t = y_{t-1}. The polynomial
% matrix D is given in condensed form in matrix DA.
%
% Input arguments:
%          y: an m x s matrix
%         DA= matrix of the form [DAr Indxr DAs Indxs], where DAr and DAs
%             are the parameterizations of the regular and seasonal
%             differencing matrix polynomials, and Indxr and Indxs are two
%             index vectors to identify the l.i. rows of DAr and DAs.
%         nr= number of regular unit roots
%         ns= number of seasonal unit roots
%       seas: seasonality
% Output arguments:
%                  yd: the series D(B)*y_t
%                  Dr: regular differencing matrix polynomial
%                  Ds: seasonal differencing matrix polynomial
%              ferror: a flag for erros
%
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
Dr = [];
Ds = [];
[ny, s] = size(y);
[nd, md] = size(DA);
yd = y;

mds = md;
if (ns > 0)
    Ds(:, :, 1) = eye(s);
    Indxs = DA(:, md);
    sncols = s - sum(Indxs);
    mdm1 = md - 1;
    mds = mdm1 - sncols;
    betaors = DA(:, mdm1-sncols+1:mdm1);
    Ds(:, :, seas+1) = -betaors * pinv(betaors'*betaors) * betaors';
    [yd, ferror] = dtimesy(Ds, yd);
end
if (nr > 0)
    Dr(:, :, 1) = eye(s);
    if s == 1
        nda = size(DA, 1);
        Dr(:, :, 2:nda) = DA(2:nda, 1);
    else
        Indxr = DA(:, mds);
        rncols = s - sum(Indxr);
        mdsm1 = mds - 1;
        betaor = DA(:, mdsm1-rncols+1:mdsm1);
        Dr(:, :, 2) = -betaor * pinv(betaor'*betaor) * betaor';
    end
    [yd, ferror] = dtimesy(Dr, yd);
end
