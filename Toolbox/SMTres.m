function [smtres] = SMTres(out)
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
% res: a vector containing the OLS residuals
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
    Y = out.model.Y;
    beta = out.model.hb;
    s = out.freq;
    d = out.model.d;
    ds = out.model.ds;
    ndelta = d + ds * s;
    %
    %compute initial conditions
    %
    [ins, ii, ferror] = incossm(T, H, ndelta);
    
    [nalpha, junk] = size(T);
    %
    % Computation with function smoothgen.m
    %
    % Function smoothgen smooths a general vector:
    % Y_t = U_t*beta + C_t*alpha_t + D_t*epsilon_t
    % In this case, it is desired to smooth:
    % Y_t = D_t*epsilon_t
    % Hence, U_t = C_t = 0 C = zeros(1,nalpha);
    U = [];
    mucd = 1;
    % C = zeros(ly*mucd,nalpha);
    % D = repmat([G;H],ly,1);
    C = zeros(1, nalpha);
    D = 1.;
    u = y - Y * beta;
    %
    %with this model, we are actually smoothing a_{t+1}. Thus, we have to
    %start at y_{0}
    %
    u = [NaN; u];
    %
    %run the disturbance smoother to compute the smoothed innovations
    %
    [DS1, MDS1] = smoothgen(u, [], Z, G, W, T, H, ins, ii, mucd, U, C, D);
    smtres = DS1(1:end-1);
else
    smtres = [];
end
