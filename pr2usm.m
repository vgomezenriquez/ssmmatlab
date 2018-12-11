function [X, Z, G, W, T, H, ins, i, ferror] = pr2usm(xx, xf, str)
% PURPOSE: given a structure containing information about a univariate
% structural model, it passes the parameters to the state space form
%---------------------------------------------------
% USAGE: [X,Z,G,W,T,H,ins,i,ferror] = pr2usm(xx,xf,str)
% where:
%          xx       = a vector containing the estimated parameters
%          xf       = a vector containing the fixed parameters
%          str      = a structure containing the initial model information
%---------------------------------------------------
%---------------------------------------------------
% RETURNS: updated matrices and initial conditions of the state space
%          form of a univariate structural model:
%
%          y_t = X_t*beta + Z_t*alpha_t + G_t*epsilon_t
%          alpha_{t+1}= W_t*beta + T_t*alpha_t + H_t*epsilon_t
%
%          where epsilon_t is (0,sigma^2I),
%
%          with initial state
%
%          alpha_1= c + W_0*beta + a_1 + A_1*delta
%
%          where c is (0,Omega) and delta is (0,kI) (diffuse)
%          More specifically:
%          X : an (n x nbeta) matrix containing the X_t matrices;
%              a  (1 x nbeta) if it is time invariant;
%              it can be []
%          Z : an (n x nalpha) matrix containing the Z_t matrices;
%              a  (1 x nalpha) matrix if it is time invariant
%          G : an (n x nepsilon) matrix containing the G_t matrices;
%              a  (1 x nepsilon) matrix if it is time invariant
%          W : an (n*nalpha x nbeta) matrix containing the W_t matrices;
%              an (nalpha x nbeta) matrix if it is time invariant;
%              it can be []
%          T : an (n*nalpha x nalpha) matrix containing the T_t matrices;
%              an (nalpha x nalpha) matrix if it time invariant
%          H : an (n*nalpha x nepsilon) matrix containing the H_t matrices;
%              an (nalpha x nepsilon) if it is time invariant
%        ins : an nalpha x (cc+cw0+ca1+cca1) matrix containing the initial
%              state information, according to array i below
%          i : a  1 x 4 array containing 4 integers, i=[cc cw0 ca1 cca1],
%              where
%              cc   = nalpha if c is not missing (0 if c missing)
%              cw0  = number of columns in W_0 (0 if W_0 missing)
%              ca1  = 1 if a_1 is not missing (0 if a_1 missing)
%              cca1 = number of columns in A_1 (0 if A_1 missing)
%     ferror : flag for errors
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

ins = str.ins;
i = str.i;
X = str.X;
W = str.W;
Z = str.Z;
s = str.freq;
pfix = str.pfix;
pvar = str.pvar;
stord = str.stord;
conc = str.conc;
trend = str.trend;
slope = str.slope;
seas = str.seas;
irreg = str.irreg;
arp = str.arp;

cycle = str.cycle;

sm1 = s - 1;
npar = length(pfix) + length(pvar);
x = zeros(1, npar);
x(pfix) = xf;
x(pvar) = xx;

xp = zeros(1, max(stord));
xp(stord) = x;
xp(conc) = 1; %concentrated parameter

%system matrices
T = [];
G = [];
H = [];
C = [];
if trend == 1
    if slope == 0
        T = 1;
        H = xp(2);
    elseif slope == 1
        T(1:2, 1:2) = [1, 1; 0, 1];
        H(1, 1) = xp(2);
        H(2, 2) = xp(3);
    elseif slope == -1
        T = 1;
        H = xp(2);
    end
elseif trend == 2
    if slope == 0
        T(1, 1) = 1;
        H(1, 1) = 2 * xp(2);
        H(2, 2) = xp(2);
        T(2, 2) = 0;
    elseif slope == 1
        den = abs(xp(3)) + abs(xp(2));
        th = (abs(xp(3)) - abs(xp(2))) / den;
        T(1:2, 1:2) = [0, 1; -1, 2];
        T(3, 3) = 0;
        H(1, 1) = (3 + th) * den;
        H(2, 1) = (5 + 3 * th) * den;
        H(3, 2) = den;
        C(3, 3) = den^2;
    elseif slope == -1
        T(1, 1) = 1;
        H(1, 1) = 2 * xp(2);
        H(2, 2) = xp(2);
        T(2, 2) = 0;
    end
end


if seas == 1 %stochastic dummy seasonality
    [nT, junk] = size(T);
    [junk, nH] = size(H);
    T(nT+1, nT+1:nT+sm1) = -ones(1, sm1);
    T(nT+2:nT+sm1, nT+1:nT+sm1-1) = eye(sm1-1);
    H(nT+1, nH+1) = xp(4);
    H(nT+sm1, nH+1) = 0;
elseif seas == 2 %trigonometric seasonality
    [nT, junk] = size(T);
    [junk, nH] = size(H);
    AA = zeros(sm1);
    for j = 1:sm1, AA(j, j) = xp(4);
    end
    T(nT+sm1, nT+sm1) = -1;
    H(nT+1:nT+sm1, nH+1:nH+sm1) = AA;
    if s == 4
        T(nT+1:nT+2, nT+1:nT+2) = [0, 1; -1, 0];
    elseif s == 12
        sqrt3d2 = sqrt(3) / 2;
        T(nT+1:nT+2, nT+1:nT+2) = [sqrt3d2, .5; -.5, sqrt3d2];
        T(nT+3:nT+4, nT+3:nT+4) = [0.5, sqrt3d2; -sqrt3d2, 0.5];
        T(nT+5:nT+6, nT+5:nT+6) = [0, 1; -1, 0];
        T(nT+7:nT+8, nT+7:nT+8) = [-0.5, sqrt3d2; -sqrt3d2, -0.5];
        T(nT+9:nT+10, nT+9:nT+10) = [-sqrt3d2, 0.5; -0.5, -sqrt3d2];
    end
elseif seas == 4 %Butterworth tan seasonality
    [nT, junk] = size(T);
    [junk, nH] = size(H);
    T(nT+sm1, nT+sm1) = -1;
    T(nT+s, nT+s) = 0;
    if s == 4
        T(nT+1:nT+2, nT+1:nT+2) = [0, 1; -1, 0];
        H(nT+2, nH+1) = -2 * xp(4);
        H(nT+3, nH+2) = -2 * xp(4);
        H(nT+4, nH+3) = sqrt(2) * xp(4);
    elseif s == 12
        sqrt3 = sqrt(3);
        T(nT+1:nT+2, nT+1:nT+2) = [0, 1; -1, sqrt3];
        T(nT+3:nT+4, nT+3:nT+4) = [0, 1; -1, 1];
        T(nT+5:nT+6, nT+5:nT+6) = [0, 1; -1, 0];
        T(nT+7:nT+8, nT+7:nT+8) = [0, 1; -1, -1];
        T(nT+9:nT+10, nT+9:nT+10) = [0, 1; -1, -sqrt3];
        H(nT+1:nT+2, nH+1) = [sqrt3 * xp(4); xp(4)];
        H(nT+3:nT+4, nH+2) = [xp(4); -xp(4)];
        H(nT+6, nH+3) = -2 * xp(4);
        H(nT+7:nT+8, nH+4) = [-xp(4); -xp(4)];
        H(nT+9:nT+10, nH+5) = [-sqrt3 * xp(4); xp(4)];
        H(nT+11, nH+6) = -2 * xp(4);
        H(nT+12, nH+7) = sqrt(6) * xp(4);
    end
    C(nT+s, nT+s) = floor(s/2) * xp(4)^2;
end

if cycle == 1 %structural model cycle (Harvey)
    [nT, junk] = size(T);
    [junk, nH] = size(H);
    xp6 = xp(6);
    sc = xp6 * sqrt((1. - xp(7)^2));
    AA = diag([sc, sc]);
    xp8 = xp(8);
    T(nT+1:nT+2, nT+1:nT+2) = xp(7) * [cos(xp8), sin(xp8); -sin(xp8), cos(xp8)];
    H(nT+1:nT+2, nH+1:nH+2) = AA;
    C(nT+1:nT+2, nT+1:nT+2) = diag([xp6^2, xp6^2]);
elseif cycle == 2 %Butterworth sine cycle
    [nT, junk] = size(T);
    [junk, nH] = size(H);
    xp6 = xp(6);
    xp7 = xp(7);
    xp8 = xp(8);
    T(nT+1:nT+2, nT+1:nT+2) = [0, 1; -xp7 * xp7, 2 * xp7 * cos(xp8)];
    H(nT+1:nT+2, nH+1) = [xp6; xp7 * cos(xp8) * xp6];
    J = H(nT+1:nT+2, nH+1);
    Tcyc = T(nT+1:nT+2, nT+1:nT+2);
    %  N=J*J'; Omega = mlyapunov(Tcyc,N);
    Omega = dlyapsq(Tcyc, J);
    Omega = Omega' * Omega;
    C(nT+1:nT+2, nT+1:nT+2) = Omega;
end


[mH, nH] = size(H);
if irreg > 0
    G(1, nH+1) = xp(1);
    H(1, nH+1) = 0;
else
    G = zeros(1, nH);
end

if arp > 0
    [nT, junk] = size(T);
    [junk, nH] = size(H);
    if arp > 1
        Tar = zeros(arp);
        Tar(1:arp-1, 2:arp) = eye(arp-1);
        Tar(arp, :) = -fliplr(xp(end-arp+1:end));
    else
        Tar = -xp(end);
    end
    T(nT+1:nT+arp, nT+1:nT+arp) = Tar;
    H(nT+arp, nH+1) = xp(5);
    G(end+1) = 0;
    J = H(nT+1:nT+arp, nH+1);
    %  N=J*J'; Omega = mlyapunov(Tar,N);
    Omega = dlyapsq(Tar, J);
    Omega = Omega' * Omega;
    C(nT+1:nT+arp, nT+1:nT+arp) = Omega;
end

[nC, mC] = size(C);
[nalpha, junk] = size(T);
if (nC > 0) & (nC < nalpha)
    C(nalpha, nalpha) = 0.;
    nC = nalpha;
end
if nC > 0
    ins(:, 1:nC) = C;
end
