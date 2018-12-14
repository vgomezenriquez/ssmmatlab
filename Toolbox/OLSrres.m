function [recrs, recr, srecr] = OLSrres(out)
%
%
% This function obtains the OLS residuals after having used function
% arimaestos, arimaestni or arimaestwi
%
% input arguments:
% out: a structure, the output of function arimaestos, arimaestni or
% arinamestwi
%
% output arguments:
%       recrs: standardized recursive residuals
%        recr: recursive residuals
%       srecr: covariance matrices of recursive residuals
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
if isfield(out.model, 'hb')
    X = out.model.matsis.X;
    W = out.model.matsis.W;
    Z = out.model.matsis.Z;
    G = out.model.matsis.G;
    H = out.model.matsis.H;
    T = out.model.matsis.T;
    lam = out.model.lam;
    y = out.orig;
    if lam == 0
        y = log(y);
    end
    s = out.freq;
    d = out.model.d;
    ds = out.model.ds;
    ndelta = d + ds * s;
    %
    %compute initial conditions
    %
    [ins, ii, ferror] = incossm(T, H, ndelta);
    
    [KKP, PT, hd, Md, initf, recrs, recr, srecr] = scakff(y, X, Z, G, W, T, H, ins, ii);
    chb = 1;
    %   [recrs2,f1,hd1,Md1,A1,P1,ne]=tskfsribf(y,X,Z,G,W,T,H,ins,ii,chb);
    %   plot(recrs-recrs2)
    %   pause
else
    recrs = [];
    recr = [];
    srecr = [];
end
