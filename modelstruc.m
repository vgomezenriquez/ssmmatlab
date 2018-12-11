function models = modelstruc(xv, y, Y, pfix, pvar, xf, modescr, datei, npr)
%************************************************************************
% This function creates a structure containing all the information
% needed for a structural model.
% This function allows for other values of freq besides 1, 4 and 12.
% Note that freq is a field of structure datei and that it is also the
% number of seasons. Therefore, it determines the seasonal component.
%
%    INPUTS:
%       xv : vector with parameters to be estimated
%           (see description below *)
%        y : data vector
%        Y : matrix for regression variables. It contains the stack
%            of the Y_t matrices
%     pfix : array with fixed parameter indices
%     pvar : array with variable parameter indices
%       xf : vector with fixed parameters (see description below *)
%  modescr : structure with the following fields:
%            (for the list of the codes of the components, see description
%            below **)
%            .trend : trend code
%            .slope : slope code
%             .seas : seasonal code
%            .cycle : cycle code
%              .xl1 : lower bound of the frequency interval in which the
%                     cycle is supposed to be defined
%              .xl2 : upper bound of the frequency interval in which the
%                     cycle is supposed to be defined
%               .ar : AR code
%            .irreg : irregular code
%            .conc : index for the parameter that is concentrated out
%                    (see description below ***)
%            .stord : array containing parameter indices
%                    (see description below ****)
%    datei : calendar structure (output of function cal)
%      npr : number of forecast
%---------------------------------------------------------------------
% * xv and xf are subvectors of the parameter vector x, where the parameters,
%   except the one that is concentrated out, are put in the following order:
%   1 - irregular standard deviation
%   2 - level standard deviation
%   3 - slope standard deviation
%   4 - seasonal standard deviation
%   5 - autoregressive standard deviation
%   6 - cycle standard deviation
%   7,8 - cycle parameters, rho and frequency
%   9,10,.. autoregressive parameters
%
% ** codes for the components:
%    trend = -1  constant
%             1  stochastic
%             2  Butterworth tangent
%    slope = -1  constant
%             1  stochastic
%    seas  = -1  fixed dummy seasonality
%             1  stochastic dummy seasonality
%             2  trigonometric seasonality
%             4  Butterworth tangent
%    cycle =  1  structural model cycle
%             2  Butterworth sine cycle
%    irreg =  1  stochastic
%      ar  =  k  autoregressive component of order k
%
% *** One of the standard deviations is concentrated out and, therefore, is not
%    estimated. The field conout contains information about this standard
%    deviation. The user can select this standard deviation or the program can
%    do it automatically instead. The biggest variance should be selected.
%
% **** stord is an index such that its i-th element indicates to which
%     component (according to the ordering above) belongs the i-th element of x.
%-------------------------------------------------------------------------
%
%   OUTPUTS:
%   models : the same structure as modescr (with .ar renamed to .arp)
%            plus the following fields:
%      matrices according to the model
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
%      More specifically:
%            .X : an (n x nbeta) matrix containing the X_t matrices;
%                 a  (1 x nbeta) if it is time invariant;
%                 it can be []
%            .Z : an (n x nalpha) matrix containing the Z_t matrices;
%                 a  (1 x nalpha) matrix if it is time invariant
%            .G : an (n x nepsilon) matrix containing the G_t matrices;
%                 a  (1 x nepsilon) matrix if it is time invariant
%            .W : an (n*nalpha x nbeta) matrix containing the W_t matrices;
%                 an (nalpha x nbeta) matrix if it is time invariant;
%                 it can be []
%            .T : an (n*nalpha x nalpha) matrix containing the T_t matrices;
%                 an (nalpha x nalpha) matrix if it time invariant
%            .H : an (n*nalpha x nepsilon) matrix containing the H_t matrices;
%                 an (nalpha x nepsilon) if it is time invariant
%     .ins : an nalpha x (cc+cw0+ca1+cca1) matrix containing the initial
%            state information, according to array i below
%       .i : a  1 x 4 array containing 4 integers, i=[cc cw0 ca1 cca1],
%            where
%            cc   = nalpha if c is not missing (0 if c missing)
%            cw0  = number of columns in W_0 (0 if W_0 missing)
%            ca1  = 1 if a_1 is not missing (0 if a_1 missing)
%            cca1 = number of columns in A_1 (0 if A_1 missing)
%
%*************************************************************************
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
%*************************************************************************

s = datei.freq;
conc = modescr.conc;
stord = modescr.stord;
if ~isfield(modescr, 'trend')
    trend = 0;
else
    trend = modescr.trend;
end
if ~isfield(modescr, 'slope')
    slope = 0;
else
    slope = modescr.slope;
end
if ~isfield(modescr, 'seas')
    seas = 0;
else
    seas = modescr.seas;
end
if ~isfield(modescr, 'cycle')
    cycle = 0;
else
    cycle = modescr.cycle;
end
if ~isfield(modescr, 'xl1')
    xl1 = 0;
else
    xl1 = modescr.xl1;
end
if ~isfield(modescr, 'xl2')
    xl2 = 0;
else
    xl2 = modescr.xl2;
end
if ~isfield(modescr, 'irreg')
    irreg = 0;
else
    irreg = modescr.irreg;
end
if ~isfield(modescr, 'ar')
    arp = 0;
else
    arp = modescr.ar;
end

ny = size(y, 1);
sm1 = s - 1;
npar = length(pfix) + length(pvar);
x = zeros(1, npar);
x(pfix) = xf;
x(pvar) = xv;

xp = zeros(1, max(stord));
xp(stord) = x;
xp(conc) = 1; %concentrated parameter

if isempty(Y)
    X = Y;
else
    X = Y(1:ny+npr, :);
end
C = [];
A = []; % matrices for initial conditions
W = [];
Z = [];
T = [];
H = [];
G = []; % system matrices

if trend == -1
    X = [ones(ny+npr, 1), X];
    if slope == -1
        X = [constant(ny+npr, 1, 1, 0, 0, 0, s), X];
    end
elseif trend == 1
    if slope == 0
        Z(1) = 1;
        T = 1;
        H = xp(2);
        A(1, 1) = 1;
    elseif slope == 1
        Z(1) = 1;
        T(1:2, 1:2) = [1, 1; 0, 1];
        H(1, 1) = xp(2);
        H(2, 2) = xp(3);
        A(1:2, 1:2) = eye(2);
        Z(2) = 0;
    elseif slope == -1
        X = [constant(ny+npr, 1, 1, 0, 0, 0, s), X];
        Z(1) = 1;
        T = 1;
        H = xp(2);
        A(1, 1) = 1;
    end
elseif trend == 2
    if slope == 0
        Z(1) = 1;
        Z(2) = 1;
        T(1, 1) = 1;
        H(1, 1) = 2 * xp(2);
        H(2, 2) = xp(2);
        T(2, 2) = 0;
        A(1, 1) = 1;
    elseif slope == 1
        den = abs(xp(3)) + abs(xp(2));
        th = (abs(xp(3)) - abs(xp(2))) / den;
        Z(1) = 1;
        Z(3) = 1;
        T(1:2, 1:2) = [0, 1; -1, 2];
        T(3, 3) = 0;
        H(1, 1) = (3 + th) * den;
        H(2, 1) = (5 + 3 * th) * den;
        H(3, 2) = den;
        A(1:2, 1:2) = eye(2);
        C(3, 3) = den^2;
    elseif slope == -1
        X = [constant(ny+npr, 1, 1, 0, 0, 0, s), X];
        Z(1) = 1;
        Z(2) = 1;
        T(1, 1) = 1;
        H(1, 1) = 2 * xp(2);
        H(2, 2) = xp(2);
        T(2, 2) = 0;
        A(1, 1) = 1;
    end
end


if seas == -1
    X = [seasdm(ny+npr, datei), X]; %fixed dummy seasonality
elseif seas == 1 %stochastic dummy seasonality
    [nT, junk] = size(T);
    [junk, nH] = size(H);
    Z(nT+1) = 1;
    [junk, nA] = size(A);
    T(nT+1, nT+1:nT+sm1) = -ones(1, sm1);
    Z(nT+sm1) = 0;
    T(nT+2:nT+sm1, nT+1:nT+sm1-1) = eye(sm1-1);
    H(nT+1, nH+1) = xp(4);
    H(nT+sm1, nH+1) = 0;
    A(nT+1:nT+sm1, nA+1:nA+sm1) = eye(sm1);
elseif seas == 2 %trigonometric seasonality
    [nT, junk] = size(T);
    [junk, nH] = size(H);
    nTT = nT;
    Z(nT+1) = 1;
    [junk, nA] = size(A);
    AA = zeros(sm1);
    for i = 1:sm1, AA(i, i) = xp(4);
    end
    H(nT+1:nT+sm1, nH+1:nH+sm1) = AA;
    per = s;
    nh = floor(per/2);
    for j = 1:nh - 1
        [nT, junk] = size(T);
        Z(nT+1) = 1;
        omega = 2 * pi * j / per;
        T(nT+1:nT+2, nT+1:nT+2) = [cos(omega), sin(omega); -sin(omega), cos(omega)];
    end
    if (mod(per, 2) == 0) && (nh == floor(per/2))
        [nT, junk] = size(T);
        Z(nT+1) = 1;
        T(nT+1, nT+1) = -1;
        Z(nT+1) = 1;
    else
        [nT, junk] = size(T);
        Z(nT+1) = 1;
        Z(nT+2) = 0;
        omega = 2 * pi * nh / per;
        T(nT+1:nT+2, nT+1:nT+2) = [cos(omega), sin(omega); -sin(omega), cos(omega)];
    end
    A(nTT+1:nTT+sm1, nA+1:nA+sm1) = eye(sm1);
elseif seas == 4 %Butterworth tan seasonality
    [nT, junk] = size(T);
    Z(nT+1) = 1;
    [junk, nH] = size(H);
    [junk, nA] = size(A);
    T(nT+sm1, nT+sm1) = -1;
    Z(nT+sm1) = 1;
    Z(nT+s) = 1;
    T(nT+s, nT+s) = 0;
    if s == 4
        T(nT+1:nT+2, nT+1:nT+2) = [0, 1; -1, 0];
        H(nT+2, nH+1) = -2 * xp(4);
        H(nT+3, nH+2) = -2 * xp(4);
        H(nT+4, nH+3) = sqrt(2) * xp(4);
    elseif s == 12
        sqrt3 = sqrt(3);
        Z(nT+3) = 1;
        Z(nT+5) = 1;
        Z(nT+7) = 1;
        Z(nT+9) = 1;
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
    A(nT+1:nT+sm1, nA+1:nA+sm1) = eye(sm1);
    A(nT+s, nA+sm1) = 0;
    C(nT+s, nT+s) = floor(s/2) * xp(4)^2;
end

if cycle == 1 %structural model cycle (Harvey)
    [nT, junk] = size(T);
    [junk, nH] = size(H);
    Z(nT+1) = 1;
    xp6 = xp(6);
    [junk, nA] = size(A);
    sc = xp6 * sqrt((1. - xp(7)^2));
    AA = diag([sc, sc]);
    xp8 = xp(8);
    T(nT+1:nT+2, nT+1:nT+2) = xp(7) * [cos(xp8), sin(xp8); -sin(xp8), cos(xp8)];
    Z(nT+2) = 0;
    H(nT+1:nT+2, nH+1:nH+2) = AA;
    C(nT+1:nT+2, nT+1:nT+2) = diag([xp6^2, xp6^2]);
    A(nT+1:nT+2, :) = zeros(2, nA);
elseif cycle == 2 %Butterworth sine cycle
    [nT, junk] = size(T);
    [junk, nH] = size(H);
    Z(nT+1) = 1;
    [junk, nA] = size(A);
    xp6 = xp(6);
    xp7 = xp(7);
    xp8 = xp(8);
    T(nT+1:nT+2, nT+1:nT+2) = [0, 1; -xp7 * xp7, 2 * xp7 * cos(xp8)];
    Z(nT+2) = 0;
    H(nT+1:nT+2, nH+1) = [xp6; xp7 * cos(xp8) * xp6];
    J = H(nT+1:nT+2, nH+1);
    Tcyc = T(nT+1:nT+2, nT+1:nT+2);
    %  N=J*J';  Omega = mlyapunov(Tcyc,N);
    Omega = dlyapsq(Tcyc, J);
    Omega = Omega' * Omega;
    C(nT+1:nT+2, nT+1:nT+2) = Omega;
    A(nT+1:nT+2, :) = zeros(2, nA);
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
    [junk, nA] = size(A);
    if arp > 1
        Tar = zeros(arp);
        Tar(1:arp-1, 2:arp) = eye(arp-1);
        Tar(arp, :) = -fliplr(xp(end-arp+1:end));
    else
        Tar = -xp(end);
    end
    T(nT+1:nT+arp, nT+1:nT+arp) = Tar;
    Z(end+arp) = 1;
    H(nT+arp, nH+1) = xp(5);
    G(end+1) = 0;
    J = H(nT+1:nT+arp, nH+1);
    %  N=J*J'; Omega = mlyapunov(Tar,N);
    Omega = dlyapsq(Tar, J);
    Omega = Omega' * Omega;
    C(nT+1:nT+arp, nT+1:nT+arp) = Omega;
    A(nT+1:nT+arp, :) = zeros(arp, nA);
end

[nC, mC] = size(C);
[nalpha, junk] = size(T);
if (nC > 0) & (nC < nalpha)
    C(nalpha, nalpha) = 0.;
    nC = nalpha;
end
[junk, nA] = size(A);
ins = [C, A];
i = [nC, 0, 0, nA];


% Z
% T
% H
% G
% A
% ins
% i
iT = find(T);
[nT, mT] = size(T);
if (iT <= .75 * nT * mT), T = sparse(T);
end
iG = find(G);
[nG, mG] = size(G);
if (iG <= .75 * nG * mG), G = sparse(G);
end
iH = find(H);
[nH, mH] = size(H);
if (iH <= .75 * nH * mH), H = sparse(H);
end

models.Z = Z;
models.T = T;
models.X = X;
models.W = W;
models.G = G;
models.H = H;
models.ins = ins;
models.i = i;
models.stord = stord;
models.conc = conc;
models.trend = trend;
models.slope = slope;
models.seas = seas;
models.cycle = cycle;
models.xl1 = xl1;
models.xl2 = xl2;
models.irreg = irreg;
models.arp = arp;
