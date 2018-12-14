function [F, e, g, M, Pevf, A, P] = frontrearfun(xx, y, xf, chb, models)
%*************************************************************************
% Auxiliary function called in usmdk3_d.m.
% Example of ''Bivariate structural time series analysis'' in the book by
% Durbin and Koopman (2201), p. 167.
%
%       This function applies to the series y the two stage Kalman filter
%       for prediction and likelihood evaluation corresponding to the model
%
%       y_t = X_t*beta + Z_t*alpha_t + G_t*epsilon_t
%       alpha_{t+1}= W_t*beta + T_t*alpha_t + H_t*epsilon_t,
%
%       where epsilon_t is (0,sigma^2I),
%
%       with initial state
%
%       alpha_1= c + W_0*beta + a_1 + A_1*delta
%
%       where c is (0,Omega) and delta is (0,kI) (diffuse). A single
%       collapse is applied to get rid of the diffuse component.
%
% INPUTS:
%        xx : array with parameters to be estimated
%         y : an (n x p) matrix of observations
%      pfix : array with fixed parameter indices
%      pvar : array with variable parameter indices
%        xf : array with fixed parameters
%       chb = 1: compute g and M
%           = 0: do not compute g and M
%    models : structure with initial model information
%
% OUTPUTS:
%        F  : residuals e mulitplied with the factor f given by scakfle2;
%             F is used for minimization of the nonlinear sum of squares
%        e  : residual vector (Q'_2*y)
%        g  : the beta estimator
%        M  : the Mse of the beta estimator
%     Pevf  : prediction error variance
%        A  : the estimated augmented state vector at the end of filtering
%        P  : the Mse of A at the end of filtering
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
%*************************************************************************

Z = models.Z;
T = models.T;
G = models.G;
H = models.H;
W = models.W;
X = models.X;
ins = models.ins;
ii = models.i;
pvar = models.pvar;
pfix = models.pfix;
stord = models.stord;
conc = models.conc; %standard deviation concentrated out

npar = length(pfix) + length(pvar);
x = zeros(1, npar);
x(pfix) = xf;
x(pvar) = xx;

xp = zeros(1, max(stord));
xp(stord) = x;
xp(conc) = 1; %concentrated parameter

[nhb, mhb] = size(H);
H(1, 1) = xp(1);
H(2, 1) = xp(2);
H(2, 2) = xp(3);
j = 5;
for i = 4:2:nhb
    H(i, j) = xp(4);
    H(i-1, j) = xp(5);
    H(i, j+1) = xp(6);
    j = j + 2;
end
G(1, mhb-1) = xp(7);
G(2, mhb-1) = xp(8);
G(2, mhb) = xp(9);


[e, f, g, M, A, P] = scakfle2(y, X, Z, G, W, T, H, ins, ii, chb);
Pevf = Z * P * Z' + G * G'; % prediction error variance (finite sample)
F = e * f;
