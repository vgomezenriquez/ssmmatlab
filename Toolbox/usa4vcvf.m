function [F, e, g, M, Xt, Pt, nmiss, initf, Pevf] = usa4vcvf(xx, y, smth, pfix, pvar, xf, chb, m)
%**************************************************************************
% Auxiliary function called in usa4vcv_d.m.
% Reference: ``Estimating Potential Output, Core Inflation
% and the NAIRU as Latent Variables'', by Rafael Domenech
% and Victor Gomez, Journal of Business and Economic Statistics (2006)
%
%        Input parameters:
%        xx   : a vector containing the varying parameters
%        y    : an (n x p) matrix of observations
%        smth : a scalar specifying whether, e.g., smoothing or filtering
%               should be performed
%        pfix : a vector containing the indices for the fixed parameters
%        pvar : a vector containing the indices for the varying parameters
%        xf   : a vector containing the fixed parameters
%        chb  : a scalar = 1, compute g and M
%                        = 0, do not compute g and M
%        m    : number of variables (0,1,2,3 or 4)
%
%        Output parameters:
%        F    : residual vector multiplied with factor f given by scakfle2;
%               used in minimization of the nonlinear sum of squares;
%               empty if smth is different from 0 and -1
%        e    : residual vector;
%               empty if smth is different from 0 and -1
%        g    : the beta estimate
%        M    : the Mse of g
%        Xt   : an (n x nalpha matrix) containing the estimated x_{t|n} if
%               smth = 1;
%               an (n x nalpha matrix) containing the estimated x_{t|t} if
%               smth is different from 0,-1 and 1;
%               empty matrix if smth = 0 or smth = -1
%        Pt   : an ((nalpha*n) x nalpha) matrix containing the
%               Mse of x_{t|n} if smth = 1;
%               an ((nalpha*n) x nalpha) matrix containing the
%               Mse of x_{t|t} if smth is different from 0,-1 and 1;
%               empty matrix if smth = 0 or smth = -1
%        nmiss: number of missing values
%        Pevf : prediction error variance (finite sample);
%               empty if smth is different from -1
%
%
% Copyright (c) 21 July 2003 by Victor Gomez
% Ministerio de Hacienda, Direccion Gral. de Presupuestos,
% Subdireccion Gral. de Analisis y P.E.,
% Alberto Alcocer 2, 1-P, D-34, 28046, Madrid, SPAIN.
% E-mail: VGomez@sepg.minhap.es
%
% The author assumes no responsibility for errors or damage resulting from the use
% of this code. Usage of this code in applications and/or alterations of it should
% be referenced. This code may be redistributed if nothing has been added or
% removed and no money is charged. Positive or negative feedback would be appreciated.
%**************************************************************************

npar = length(pfix) + length(pvar);
x = zeros(1, npar);

if m > 0
    xx = usautrcv(xx, m, pfix, pvar, xf); %parameters are untransformed
end

x(pfix) = xf;
x(pvar) = xx;
[n, p] = size(y);
my = zeros(n-4, p);
my(:, 1) = y(5:n, 1);
my(:, 2) = y(5:n, 3) - x(8) * y(4:n-1, 3);
my(:, 3) = y(5:n, 4) - x(9) * y(4:n-1, 4);
my(:, 4) = y(5:n, 2) - x(10) * y(4:n-1, 2) - x(11) * y(3:n-2, 2) - x(12) * y(2:n-3, 2) ...
    -x(13) * y(1:n-4, 2);
nmiss = sum(sum(isnan(my))); % number of missing values

nalpha = 7; %state vector dimension

Z = zeros(p, nalpha);
Z(1, 1) = 1;
Z(1, 7) = 1;
Z(2, 2) = 1 - x(8);
Z(2, 5:7) = x(16:-1:14);
Z(3, 3) = 1 - x(9);
Z(3, 6:7) = x(18:-1:17);
Z(4, 4) = 1 - sum(x(10:13));
Z(4, 7) = x(19);
%  Z
Z = sparse(Z);

T = zeros(nalpha, nalpha);
HH = zeros(nalpha, 8);
G = zeros(p, 8);

T(1:p, 1:p) = eye(p);
T(p+1:p+2, p+2:nalpha) = eye(2);
T(nalpha, 6) = -x(20)^2;
T(nalpha, nalpha) = 2 * x(20) * cos(x(21));

HH(1:p, 1:p) = diag(x(1:4));
HH(nalpha, 8) = x(5);
H = zeros((n - 4)*nalpha, 8);
for i = 1:n - 4
    ip = (i - 1) * nalpha + 1:i * nalpha;
    HW = HH;
    if i >= 101 %105  = 1972:I  (4 data points are lost)
        HW(4, 4) = x(23);
    end
    if i >= 146 %150 = 1983:II  (4 data points are lost)
        HW(nalpha, 8) = x(22);
        HW(4, 4) = x(24);
    end
    H(ip, :) = HW;
end


% G(2:p,5:7)=diag(x(5:7));
G(2:3, 5:6) = diag(x(6:7));
G(4, 7) = 1;


%  T
%  G
%  H

T = sparse(T);
G = sparse(G);
H = sparse(H);


% X=[];
nbeta = 4;
X = zeros((n - 4)*p, nbeta);
X(17*p, 2) = 1; %21 = 1951:I   (4 data points are lost)
X(8*p, 3) = 1; %12 = 1948:IV  (4 data points are lost)
for i = 1:30
    X((8 + i)*p, 3) = .7^i; %15 = 1949:III (4 data points are lost)
end
X(19*p, 4) = 1; %23 = 1951:III (4 data points are lost)
% X(11*p,5)=1;  %15 = 1949:III (4 data points are lost)
% for i=1:30
%  X((11+i)*p,5)=.7^i;  %15 = 1949:III (4 data points are lost)
% end
% X
X = sparse(X);
W = zeros(nalpha, nbeta);
W(1, 1) = 1;
%W
W = sparse(W);


W0 = zeros(nalpha, nbeta);
W0(1, 1) = 1;
%W0
% W0=sparse(W0);

A = [eye(p); zeros(nalpha-p, p)];
% A=sparse(A);
F = T(p+1:nalpha, p+1:nalpha);
J = H(p+1:nalpha, 8);
N = J * J';
Omega = mlyapunov(full(F), full(N));
bOmega = zeros(nalpha, nalpha);
bOmega(p+1:nalpha, p+1:nalpha) = Omega;

ins = [bOmega, W0, A];
i = [nalpha, nbeta, 0, p];

% % big k specification
% bOmega(1:p,1:p)=1e+8*eye(p);
% ins=[bOmega W0];
% i=[nalpha nbeta 0 0 ];


initf = 0;
Pevf = [];
if smth == 0
    [e, f, g, M] = scakfle2(my, X, Z, G, W, T, H, ins, i, chb);
    F = e * f;
    Xt = [];
    Pt = [];
elseif smth == -1
    [e, f, g, M, A, P] = scakfle2(my, X, Z, G, W, T, H, ins, i, chb);
    F = e * f;
    Xt = [];
    Pt = [];
    Pevf = Z * P * Z' + G * G'; % prediction error variance (finite sample)
elseif smth == 1
    [Xt, Pt, g, M] = scakfs(my, X, Z, G, W, T, H, ins, i);
    F = [];
    e = [];
else
    [Xt, Pt, g, M, initf] = scakff(my, X, Z, G, W, T, H, ins, i); % regression parameters recursively estimated
    [Xt, Pt] = scakfff(my, X, Z, G, W, T, H, ins, i, g); % regression parameters fixed
    F = [];
    e = [];
end
