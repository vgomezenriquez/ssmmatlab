function [np, X, Z, G, W, T, H, ins, ii] = susmspbp(compbp, x, pfix, pvar, xf, stordt, conc, n, r, l)
%**************************************************************************
% Function to set up the state space model corresponding to the structural
% model
%
%      y_t = mu_t + epsilon_t,
%
% where
%
%      mu_{t+1} = mu_{t} + beta_{t} + zeta_{t}
%    beta_{t+1} = beta_{t} + eta_{t},
%
% plus a band-pass filter, applied as described in "The Use of Butterworth
% Filters for Trend and Cycle Estimation in Economic Time Series", Gómez,
% V. (2001), Journal of Business and Economic Statistics, 19, 365-373.
% This model is used in the paper "Monthly US Business Cycle Indicators: A
% new Multivatiate Approach Based on a Band-pass Filter", Marczak and
% Gómez, Empirical Economics, 52, 1379-1408.
%
%  INPUTS:

%         x : array with estimated parameters
%      pfix : array with fixed parameter indices
%      pvar : array with variable parameter indices
%        xf : array with fixed parameters
%    stordt : array with indices for standard deviations
%      conc : index for the standard deviation to be concentrated out
%         n : number of variables in the data matrix y
%         r : number of slope factors
%         l : numbef of quarterly series
%
%    OUTPUTS:
%        np : number of rows of matrix Tp (see below)
%        the following arguments are for generated the state space model
%        X    : an (n*p x nbeta) matrix containing the X_t matrices;
%               a  (p x nbeta) if it is time invariant;
%               it can be []
%        Z    : an (n*p x nalpha) matrix containing the Z_t matrices;
%               a  (p x nalpha) matrix if it is time invariant
%        G    : an (n*p x nepsilon) matrix containing the G_t matrices;
%               a  (p x nepsilon) matrix if it is time invariant
%        W    : an (n*nalpha x nbeta) matrix containing the W_t matrices;
%               an (nalpha x nbeta) matrix if it is time invariant;
%               it can be []
%        T    : an (n*nalpha x nalpha) matrix containing the T_t matrices;
%               an (nalpha x nalpha) matrix if it time invariant
%        H    : an (n*nalpha x nepsilon) matrix containing the H_t matrices;
%               an (nalpha x nepsilon) if it is time invariant
%        ins: an nalpha x (cc+cw0+ca1+cca1) matrix containing the initial
%             state information, according to array i below
%        ii   : a  1 x 4 array containing 4 integers, i=[cc cw0 ca1 cca1],
%             where
%             cc   = nalpha if c is not missing (0 if c missing)
%             cw0  = number of columns in W_0 (0 if W_0 missing)
%             ca1  = 1 if a_1 is not missing (0 if a_1 missing)
%             cca1 = number of columns in A_1 (0 if A_1 missing)
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
%*************************************************************************
modelsm = modstr_mulcycuswcv(x, pfix, pvar, xf, stordt, conc, n, r, l);
[junk, nthetaz] = size(compbp.den);
thetaz = zeros(n, n, nthetaz);
for i = nthetaz:-1:1
    A = repmat(compbp.den(i), n, 1);
    thetaz(:, :, nthetaz-i+1) = diag(A);
end
[junk, ndeltaz] = size(compbp.Alpha);
deltaz = zeros(n, n, ndeltaz);
for i = ndeltaz:-1:1
    A = repmat(compbp.Alpha(i), n, 1);
    deltaz(:, :, ndeltaz-i+1) = diag(A);
end
Di = compbp.Di;
num = compbp.num;
for i = 1:2
    num = deconv(num, [-1., 1.]);
end
[junk, nnumc] = size(num); %num = (1-z)^(Di-2)(1+z)^Di
numc = zeros(n, n, nnumc);
for i = nnumc:-1:1
    A = repmat(num(i), n, 1);
    numc(:, :, nnumc-i+1) = diag(A);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%smooth trend model: [(I-B*I)^2]*thetaz(B) p_t = deltaz(B)*(I+Theta*B)P_t
%
%The model for p_t is implemented in cascade form:
%     p_t = [thetaz(B)^{-1}deltaz(B)][(I-B*I)^{-2}(I+Theta*B)] P_t
%
%model for z_t based on the original form
%Model for the trend in original form.
Gz = modelsm.Gm; %Gm is zero and needs no multiplication by Lambda/sa
Tz = modelsm.Tm;
Hz = modelsm.Hm .* (compbp.Lambda / compbp.sa);
Zz = modelsm.Zm;
%model for  p_t = [thetaz(B)^{-1}deltaz(B)] z_t
[nt, mt, dp1t] = size(thetaz);
dt = dp1t - 1;
dtbn = dt * n;
dp1tbn = dp1t * n;
Ty = zeros(dp1tbn);
Ty(1:dtbn, n+1:dp1tbn) = eye(dtbn);
Hy = zeros(n*dp1t, n);
[J, ierror] = ptransfer(thetaz, deltaz, ndeltaz);
bphi = zeros(n, n, dp1t+1);
bphi(:, :, 1:dp1t) = thetaz;
for i = 1:dp1t
    Hy((i - 1)*n+1:i*n, :) = J(:, :, i);
    Ty(end-n+1:end, (i - 1)*n+1:i*n) = -bphi(:, :, dp1t-i+2);
end
Zy = zeros(n, dp1tbn);
Zy(1:n, 1:n) = eye(n);
%  The following state space form is a cascade implementation for p_t
%
%        x_{t}   = [ T_y   H_y*Z_z]x_{t-1}  + [ H_y*G_z ]e_t
%                  [  0    T_z    ]           [ H_z ]
%          y_t   = [ Z_y     0    ]x_t ,
%
[nTz, mTz] = size(Tz);
Tp = [Ty, Hy * Zz; zeros(nTz, dp1tbn), Tz];
Hp = [Hy * Gz; Hz];
Zp = [Zy, zeros(n, mTz)];
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%model for c_t = thetaz(B)^{-1}[numc(B)(I+Theta*B)] C_t based on the
%original form
K = modelsm.Tm(1:n, end);
Tw = zeros(n);
Zw = eye(n);
Dp5zeta = modelsm.Hm(1:n, 1:n) .* (1 ./ compbp.sa);
Dp5eta = modelsm.Hm(n+1:n+r, n+1:n+r) .* (1 ./ compbp.sa);
Hw = [-Dp5zeta, zeros(n, r)];
Gw = [Dp5zeta, K * Dp5eta];
th1(:, :, 1) = Gw;
th1(:, :, 2) = Hw;
[th, ierror] = pmatmul(numc, th1);
[Zc, Tc, Hc, ierror] = qarmax2ss12(thetaz, th);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%monthly model including smooth trend and cycle
% matrix Tmpc
[np, junk] = size(Tp);
[nc, junk] = size(Tc);
nalpha = np + nc;
Tmpc = zeros(nalpha);
Tmpc(1:np, 1:np) = Tp;
Tmpc(np+1:np+nc, np+1:np+nc) = Tc;
% matrix Hmpc
[nhp, mhp] = size(Hp);
[nhc, mhc] = size(Hc);
mhppmhc = mhp + mhc;
mh = mhppmhc + n;
Hmpc = zeros(nalpha, mh);
Hmpc(1:nhp, 1:mhp) = Hp;
Hmpc(nhp+1:nhp+nhc, mhp+1:mhppmhc) = Hc;
% matrix Gmpc
np1 = n + 1;
tn = 2 * n;
stord = setdiff(stordt, conc);
xs = zeros(1, length(stordt));
xs(stord) = x(n:end);
xs(conc) = 1.;
Gmpc = zeros(n, mh);
Gmpc(:, mhppmhc+1:end) = diag(xs(np1:tn)');
% matrices Zmpc, Wmpc and Xmpc
Zmpc = zeros(n, nalpha);
Zmpc(1:n, 1:n) = eye(n);
Zmpc(1:n, np+1:np+n) = eye(n);
Wmpc = [];
Xmpc = [];

T = Tmpc;
H = Hmpc;
Z = Zmpc;
G = Gmpc;

%initial conditions
%zeta_1
zeta1 = dlyapsq(Tc, Hc);
zeta1 = zeta1' * zeta1;
%eta_1
[ndelta, junk] = size(Tz);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[insp, ip, ferror] = incossm(Tp, Hp, ndelta);
A = insp(:, np+1:end);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Us = Tp(1:dp1tbn,1:dp1tbn);
% Un = Tp(dp1tbn+1:end,dp1tbn+1:end);
% U12 = Tp(1:dp1tbn,dp1tbn+1:end);
% [XX,ferror] = mclyapunov(Us,Un,U12);
% A=[-XX; eye(ndelta)];
% Hs=Hp(1:dp1tbn,:) + XX*Hp(dp1tbn+1:end,:);
% V=dlyapsq(Us,Hs); V=V'*V;
% Rs=[eye(dp1tbn); zeros(ndelta,dp1tbn)]; V=Rs*V*Rs';
% [nalpha,junk]=size(Tp);
% ip=[nalpha 0 0 ndelta]; insp=[V A];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sp = 1:ip(1);
p = insp(:, sp);
% alpha_1 = A*delta + c
nalpha = np + nc;
Aa = zeros(nalpha, ndelta);
Aa(1:np, :) = A;
c = zeros(nalpha);
c(1:np, 1:np) = p;
c(np+1:end, np+1:end) = zeta1;
ii = [nalpha, 0, 0, ndelta];
ins = [c, Aa];
X = [];
W = [];
