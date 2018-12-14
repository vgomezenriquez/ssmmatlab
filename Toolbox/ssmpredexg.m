function [ypr, mypr, alpr, malpr] = ssmpredexg(n, x, stx, sts)
%
%        This function computes n forecasts of the exogenous part of a
%        VARMAX model in echelon form.
%
% The state space echelon form is:
%
%  alpha_{t+1} = F*alpha_{t} + B*x_t{t} + K*a_{t}
%      y_{t}   = Y_{t}*beta + H*alpha_{t} + D*x_{t}  + a_{t}
%
% Thus,
%
%      y_{t} = Y_{t}*beta + V_{t} + U_{t},
%
% where V_{t} is the exogenous part,
%
%      V_{t} = [H*(zI - F)_{t}^{-1}B*x_{t} + D*x_{t}] + H*F^{t-1}*m1
%
%---------------------------------------------------
% USAGE: [ypr,mypr,alpr,malpr]=ssmpredexg(n,x,stx,sts)
%
%        Input parameters:
%        n  : number of forecasts
%        x  : a (nobs x mx) matrix containing the input variables
%      stx  : a structure with fields .T, .Z, .B and .D such that
%             T=F, Z=H, and B and D are matrices of the VARMAX model
%      stx  = a structure. If it is empty, the inputs are nonstochastic and
%             their forecasts are included in x. If it is not empty, it has
%             fields .T, .Z, .H and .G corresponding to an VARMA model
%             followed by the inputs.
%
%        Output parameters:
%        ypr : a (p x n) matrix containing the forecasts of the
%              observations
%        mypr: a (p x p x n) array containing the Mse of forecasts
%        alpr: an (nalpha x n) matrix containing the forecasts of the state
%              vector
%       malpr: an (nalpha x nalpha x n) array containing the Mse of alpr
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

[nobs, junk] = size(x);
%compute the exogenous part
F = stx.T;
H = stx.Z;
B = stx.B;
D = stx.D;
m1 = stx.m1;
[neqs, nalpha] = size(H);
inc = 1;
krov = stx.kro;
[z, rx1] = varmafil(x, F, H, B, D, krov, inc);
yfvec = vec(z'); %put the observations in a vector
yfvec = yfvec + rx1 * m1; %add estimated initial state
yx = reshape(yfvec, neqs, nobs);
V = yx'; %transform from vector to matrix

if isempty(sts)
    %inputs are deterministic. Forecasts are included in x.
    ypr = V(end-n+1:end, :)';
    mypr = [];
    alpr = [];
    malpr = [];
else
    %inputs are stochastic. The input model is included in structure sts.
    Fz = sts.T;
    Hz = sts.Z;
    Kz = sts.H;
    Jz = sts.G;
    W = [];
    X = [];
    [mx, ngamma] = size(Hz);
    %set up combined model
    T = [F, B * Hz; zeros(ngamma, nalpha), Fz];
    HH = [B * Jz; Kz];
    Z = [H, D * Hz];
    G = D * Jz;
    %initial conditions
    g1 = [m1; zeros(ngamma, 1)];
    Sigma = dlyapsq(T, HH);
    Sigma = Sigma' * Sigma;
    ins = [Sigma, g1];
    i = [nalpha + ngamma, 0, 1, 0];
    %run Kalman filter
    chb = 1;
    [e, f, hb, Mb, A, P] = scakfle2(V, X, Z, G, W, T, HH, ins, i, chb);
    %compute forecasts
    [ypr, mypr, alpr, malpr] = ssmpred(n, neqs, A, P, X, Z, G, W, T, HH, hb, Mb);
end
