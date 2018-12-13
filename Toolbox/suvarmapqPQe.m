function [str, ferror] = suvarmapqPQe(Lambda, alpha, betap, th, Th, Sigma, freq)
%**********************************************************************
% PURPOSE: sets up a VARMA model given the matrix polynomials and the
% frequency (number of observations per year). It returns a structure
% containing the model information
%---------------------------------------------------
% USAGE: [str,ferror] = suvarmaxpqPQ(phi,th,Phi,Th,Sigma,freq)
% where:   Lambda   = the regular AR matrix polynomial in ecf
%          alpha    = a matrix such that Pi = alpha*betap
%          betap    = a matrix such that Pi = alpha*betap
%          th       = the regular MA matrix polynomial
%          Th       = the seasonal MA matrix polynomial
%          Sigma    = the innovations covariance matrix
%          freq     = number of observations per year
%---------------------------------------------------
% RETURNS: str = a structure containing the following fields
%            .s: dimension of Y_t
%         .freq: number of observations per year
%          .phi: phi(z) matrix polynomial
%           .th: theta(z) matrix polynomial
%          .Phi: Phi(z^freq) matrix polynomial
%           .Th: Theta(z^freq) matrix polynomial
%        .Sigma: innovations covariance matrix, where the element
%               (1,1) has been concentrated out of the likelihood
%           .Lh: vech of the Cholesky factor of Sigma
%         .phin: phi(z) matrix polynomial, where each parameter to estimate
%                is replaced with NaN and each fixed parameter is replaced
%                with zero
%         .Phin: Phi(z^freq) matrix polynomial, where each parameter to estimate
%                is replaced with NaN and each fixed parameter is replaced
%                with zero
%          .thn: Theta(z^freq) matrix polynomial, where each parameter to estimate
%                is replaced with NaN and each fixed parameter is replaced
%                with zero
%          .Thn: phi(z) matrix polynomial, where each parameter to estimate
%                is replaced with NaN and each fixed parameter is replaced
%                with zero
%          .Lhn: Lh vector, where each parameter to estimate is replaced with NaN
%                and each fixed parameter is replaced with zero
%            .x: parameter vector, including fixed and variable
%                parameters
%        .nparm: number of parameters
%           .xi: array of ones and zeros, where each one indicates a parameter in
%                x that is to be estimated and each zero indicates a fixed parameter
%           .xv: array of parameters to estimate
%           .xf: array of fixed parameters
%        .phirs: product of phi(z) and Phi(z^freq)
%         .thrs: product of theta(z) and Theta(z^freq)
%         .nr  : number of unit roots in the model
%         .ns  : number of seasonal unit roots in the model (not used)
%         .xd  : parameter vector for betaor (the unit root part), including
%                fixed and variable parameters
%       .nparmd: number of parameters for the unit root part
%         .xid : array of ones and zeros, as in xi, for the unit root part
%         .xvd : array of parameters to estimate for the unit root part
%         .xfd : array of fixed parameters for the unit root part
%
%      Z,T,H,G,X,W matrices of the state space model
%
%                 Y_t = X_t*beta + Z_t*alpha_t + G_t*epsilon_t
%          alpha_{t+1}= W_t*beta + T_t*alpha_t + H_t*epsilon_t
%
%          where epsilon_t is (0,sigma^2I),
%
%          with initial state
%
%          alpha_1= c + W_0*beta + a_1 + A_1*delta
%
%          where c is (0,Omega) and delta is (0,kI) (diffuse)
%
%      More specifically:
%            .Z: a (s x nalpha) matrix
%            .T: a (nalpha x nalpha) matrix
%            .H: a (nalpha x nepsilon) matrix
%            .G: a (s x nepsilon) matrix
%            .X: []
%            .W: []
%
%     Note: user can subsequently incorporate regression variables into the
%           state space model given by suvarmapqPQe by an appropriate
%           specification of matrix X
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
%*************************************************************************

ferror = 0;
str = [];

[np, mp, pm1] = size(Lambda);
[nP, ma] = size(alpha);
[nb, mP] = size(betap);
[nt, mt, qm1] = size(th);
[nT, mT, Qm1] = size(Th);
[ns, ms] = size(Sigma);

if ((np ~= nt) || (np ~= nP) || (np ~= nT) || (np ~= ns) || ...
        (ma ~= nb) || ...
        (mp ~= mt) || (mp ~= mP) || (np ~= mT) || (np ~= ms))
    ferror = 1;
    disp('dimensions of Lambda, alpha, betap, th, Th and Sigma')
    disp('should agree in suvarmaxpqPQe');
    return
end

if nargin < 7
    ferror = 2;
    disp('the number of arguments in suvarmaxpqPQe is seven');
    return
end

str.s = np;
str.freq = freq;
str.Lambda = Lambda;
str.alpha = alpha;
str.betap = betap;
str.th = th;
str.Th = Th;
str.Sigma = Sigma;


[R, pp] = chol(Sigma);
L = R';
if (pp > 0)
    disp('covariance matrix of residuals singular in suvarmaxpqPQ')
    ferror = 3;
    return
end
L = L ./ L(1, 1); %sigma(1,1) is concentrated out of the
%likelihood
Lh = vech(L);
str.Lh = Lh;

Lambdan = Lambda;
Lambdan(:, :, 2:end) = NaN(size(Lambda(:, :, 2:end)));
thn = th;
thn(:, :, 2:end) = NaN(size(th(:, :, 2:end)));
alphan = NaN(size(alpha));
Thn = Th;
Thn(:, :, 2:end) = NaN(size(Th(:, :, 2:end)));
str.Lambdan = Lambdan;
str.alphan = alphan;
str.thn = thn;
str.Thn = Thn;
Lhn = Lh;
Lhn(2:end) = NaN(size(Lh(2:end)));
str.Lhn = Lhn;

%create vector containing model parameters
x = [];
for i = 2:pm1
    x = [x, vec(Lambda(:, :, i))'];
end
x = [x, vec(alpha)'];
for i = 2:qm1
    x = [x, vec(th(:, :, i))'];
end
for i = 2:Qm1
    x = [x, vec(Th(:, :, i))'];
end
x = [x, Lh(2:end)'];
nparm = length(x);
%  x
str.x = x;
str.nparm = nparm;
xi = ones(size(x));
xv = x;
xf = [];
str.xi = xi;
str.xv = xv;
str.xf = xf;

%transform error correction model to model in differences
[phi, D, DA, ferror] = mecf2mid(Lambda, alpha, betap);
str.phi = phi;
str.D = D;
str.DA = DA;
Phi = phi(:, :, 1);
str.Phi = Phi;
[phirs, thrs, H, F, G, J, ferror] = varmapqPQ2ssm(phi, th, Phi, Th, L, str);
str.phirs = phirs;
str.thrs = thrs;
str.Z = H;
str.T = F;
str.H = G;
str.G = J;
str.X = [];
str.W = [];
%add unit root information to the model structure (parameters corresponding
%to betapor)
nr = np - ma; %number of unit roots in the model
[str, ferror] = aurirvarmapqPQ(str, nr, DA);
