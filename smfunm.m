function [F, e, g, M, Pevf, A, P] = smfunm(xx, y, pfix, pvar, xf, chb, models)
%*************************************************************************
%  This function evaluates the residuals for the nonlinear minimization
%  of the sum of squares of a canonical structural model
%
%   INPUTS:
%      xx : vector with parameters to be estimated
%       y : data vector
%    pfix : array with fixed parameter indices
%    pvar : array with variable parameter indices
%      xf : vector with fixed parameters
%     chb = 1 : compute beta estimator hb and MSE of hb
%         = 0 : do not compute hb and MSE of hb
%  models : structure containing model information
%
%  OUTPUTS:
%       F : residual vector multiplied with factor f given by scakfle2;
%           used in minimization of nonlinear sum of squares
%       e : residual vector
%       g : the beta estimator
%       M : the MSE of the beta estimator
%    Pevf : prediction error variance
%       A  : the estimated augmented state vector at the end of filtering
%       P  : the MSE of A at the end of filtering
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

models.pfix = pfix;
models.pvar = pvar;
[X, Z, G, W, T, H, ins, i, ferror] = pr2usmm(xx, xf, models);

if isfield(models, 'sqrtfil')
    sqrtfil = models.sqrtfil;
else
    sqrtfil = 0;
end
if (sqrtfil == 1)
    %use the following function instead of scakfle2 in case of numerical problems
    [e, f, g, M, A, LP] = scakflesqrt(y, X, Z, G, W, T, H, ins, i, chb);
    P = LP * LP';
else
    [e, f, g, M, A, P] = scakfle2(y, X, Z, G, W, T, H, ins, i, chb);
end
Pevf = Z * P * Z' + G * G'; % prediction error variance (finite sample)
F = e * f;
