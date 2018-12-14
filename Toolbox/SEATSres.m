function [res, str] = SEATSres(out)
%
%
% This function obtains the residuals given by program SEATS by Gomez and
% Maravall (see Gomez, V., and Maravall, A., 2001. Programs TRAMO and
% SEATS, Instructions for the User (Beta Version: June 1997) (Working Paper
% No. 97001). Direccion General De Presupuestos, Ministry of Finance,
% Madrid, Spain.)
%
% input arguments:
% out: a structure, the output of function arimaestos
%
% output arguments:
% res: a vector containing the SEATS residuals
% str: a structure containing the pure MA model used in Burman (1980)
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
y = out.orig;
s = out.freq;
if isfield(out.model, 'nreg')
    beta = out.model.hb;
    Y = out.model.Y;
else
    Y = [];
end

lam = out.model.lam;
if lam == 0
    y = log(y);
end
yd = y;
Yd = Y;
if out.model.d > 0
    d = out.model.d;
    for i = 1:d
        yd = diferm(yd, 1);
        if ~isempty(Yd)
            Yd = diferm(Yd, 1);
        end
    end
end
if out.model.ds > 0
    ds = out.model.ds;
    for i = 1:ds
        yd = diferm(yd, s);
        if ~isempty(Yd)
            Yd = diferm(Yd, s);
        end
    end
end

if ~isempty(Yd)
    yd = yd - Yd * beta;
end

if out.model.p > 0
    p = out.model.p;
    phi = out.model.phi;
    phi = [1., phi];
else
    p = 0;
    phi = 1.;
end
if out.model.ps > 0
    ps = out.model.ps;
    phis = out.model.phis;
    phis = [1., zeros(1, s-1), phis]; %deg(phis) <= 1
else
    ps = 0;
    phis = 1.;
end
pt = p + ps * s;
if pt > 0
    phic = fliplr(phi);
    phisc = fliplr(phis);
    phitc = conv(phic, phisc);
    phit = fliplr(phitc);
    phit = sparse(phit);
    [n, m] = size(yd);
    u = zeros(n-pt, 1);
    for i = pt + 1:n
        u(i-pt) = yd(i) + phit(2:end) * yd(i-1:-1:i-pt);
    end
else
    u = yd;
end
%
%now we compute the smoothed innovations using the Kalman filter and the
%disturbance smoother with the pure MA model given by Burman (1980).
%
if out.model.q > 0
    q = out.model.q;
    th = out.model.th;
    th = [1., th];
else
    q = 0;
    th = 1.;
end
if out.model.qs > 0
    qs = out.model.qs;
    ths = out.model.ths;
    ths = [1., ths]; %deg(ths) <= 1
else
    qs = 0;
    ths = 1.;
end
if q + qs > 0
    %
    % put pure MA model into state space form
    %
    Sigmaf = out.model.resinf.conp;
    phir(:, :, 1) = 1.;
    Phi(:, :, 1) = 1.;
    theta(:, :, 1) = 1;
    thetas(:, :, 1) = 1;
    for i = 2:q + 1
        theta(:, :, i) = th(i);
    end
    for i = 2:qs + 1
        thetas(:, :, i) = ths(i);
    end
    [str, ferror] = suvarmapqPQ(phir, theta, Phi, thetas, Sigmaf, s);
    %
    % Computation with function smoothgen.m
    %
    % Function smoothgen smooths a general vector:
    % Y_t = U_t*beta + C_t*alpha_t + D_t*epsilon_t
    % In this case, it is desired to smooth:
    % Y_t = D_t*epsilon_t
    % Hence, U_t = C_t = 0
    X = str.X;
    W = str.W;
    Z = str.Z;
    G = str.G;
    H = str.H;
    T = str.T;
    ndelta = 0;
    %
    %compute initial conditions
    %
    [ins, ii, ferror] = incossm(T, H, ndelta);
    
    [nalpha, junk] = size(T);
    U = [];
    mucd = 1;
    % C = zeros(ly*mucd,nalpha);
    % D = repmat([G;H],ly,1);
    C = zeros(1, nalpha);
    D = 1.;
    %
    %make initial values missing to comform with Burman (1980)
    %
    ux = [NaN(q+qs*s, 1); u];
    %
    %run the disturbance smoother to compute the smoothed innovations
    %
    [DS1, MDS1] = smoothgen(ux, X, Z, G, W, T, H, ins, ii, mucd, U, C, D);
    res = DS1;
else
    res = u;
end
