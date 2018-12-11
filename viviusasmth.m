function [Xt, Pt, g, M] = viviusasmth(xx, y, pfix, pvar, xf, stordt, conc, n, r)
%*************************************************************************
% Auxiliary function called in mulcycus.m.
%
%       This function applies to the series y the two stage Kalman filter
%       for smoothing corresponding to the model
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
%         y : an (n x p) matrix of observations
%        xx : array with parameters to be estimated
%      pfix : array with fixed parameter indices
%      pvar : array with variable parameter indices
%        xf : array with fixed parameters
%       chb = 1: compute g and M
%           = 0: do not compute g and M
%    models : structure with initial model information
%         Q : orthogonal matrix needed in modstr_agtrimanbs.m for
%             tranformation of the model for stacked observations;
%         p : number of variables in the data matrix y
%        nt : length of y
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

models = modstr_viviusa(y, xx, pfix, pvar, xf, stordt, conc, n, r);
Z = models.Z;
T = models.T;
G = models.G;
H = models.H;
W = models.W;
X = models.X;
ins = models.ins;
i = models.i;

% [Xt,Pt,g,M]=scakfs(y,X,Z,G,W,T,H,ins,i);
[Xt, Pt, g, M] = scakfssqrt(y, X, Z, G, W, T, H, ins, i);
