function [Xt, Pt, g, M] = agtrimanbssmt(xx, y, pfix, pvar, xf, models, Q, p, ny)
%*************************************************************************
% Auxiliary function called in agtrimanssbs_d.m.
% Reference: ''A new State Space Methodology to
% Disaggregate Multivariate Time Series'', Journal of Time Series Analysis,
% 30, 97-124, by Gomez and Aparicio, (2009)
%
%        This function applies to the series y the augmented Kalman
%        filter and smoother corresponding to the model
%
%        y_t = X_t*beta + Z_t*alpha_t + G_t*epsilon_t
%        alpha_{t+1}= W_t*beta + T_t*alpha_t + H_t*epsilon_t
%
%        where epsilon_t is (0,sigma^2I),
%
%        with initial state
%
%        alpha_1= c + W_0*beta + a_1 + A_1*delta
%
%        where c is (0,Omega) and delta is (0,kI) (diffuse). A single
%        collapse is applied to get rid of the diffuse component.
%
% INPUTS:
%        xx : array with parameters to be estimated
%         y : an (n x p) matrix of observations
%      pfix : array with fixed parameter indices
%      pvar : array with variable parameter indices
%        xf : array with fixed parameters
%    models : structure with initial model information
%         Q : orthogonal matrix needed in modstr_agtrimanbs.m for
%             tranformation of the model for stacked observations;
%         p : number of variables in the data matrix y
%        ny : length of y
%
% OUTPUTS:
%        Xt : an (n x nalpha) matrix containing the smoothed state vector
%        Pt : an (n*nalpha x nalpha) matrix containing the Mse of Xt
%        g  : the beta estimator
%        M  : the Mse of the beta estimator
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

stcs = 0;
[modelsn, H_p, J_p, Q] = modstr_agtrimanbs(xx, pfix, pvar, xf, models, Q, stcs, p, ny);
Z = modelsn.Z;
T = modelsn.T;
G = modelsn.G;
H = modelsn.H;
W = modelsn.W;
X = modelsn.X;
ins = modelsn.ins;
i = modelsn.i;

[Xt, Pt, g, M] = scakfs(y, X, Z, G, W, T, H, ins, i);
