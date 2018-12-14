function [Z, T, H, A, Sigma, Xi] = arimam(phi, alpha, th)
%
% this function sets up the state space matrices necessary for the
% Kalman filter given the ARIMA polynomials. The state space model is
%
%    x_{t} = T x_{t-1} + H a_{t}
%    y_{t} = Z x_{t},
%
% where the initial state vector is
%
%        x_{d+1} = A*\delta + \Xi*c,
%
% and Var(c) = Sigma. See Gomez and Maravall (1994), "Estimation,
% Prediction and Interpolation for Nonstationary Series with the
% Kalman Filter", Journal of the American Statistical Association,
% 89, 611-624. The filter is initialized at time t = d+1, where d is
% the differencing degree, and the first d observations are stacked
% to form the \delta vector.

%
%  INPUTS:
%     phi   : an array containing the AR polynomial
%     alpha : an array containing the differencing polynomial
%     th    : an array containing the MA polynomial
%
%  OUTPUTS:
%        Z: the Z matrix
%        T: the T matrix
%        H: the H matrix
%        A: the A matrix
%    Sigma: the Sigma matrix
%       Xi: the Xi matrix
%
% Copyright (c) 21 July 2015 by Victor Gomez
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
%
p = length(phi) - 1;
q = length(th) - 1;
d = length(alpha) - 1;
phialph = conv(phi, alpha);
%
% dimensions for the model
%
pd = p + d;
r = max(pd, q+1);
%
% Setup matrices for the state space form
%
T = zeros(r);
T(1:r-1, 2:r) = eye(r-1);
if pd > 0
    T(r, r-p-d+1:r) = -phialph(1:pd);
end
Z = zeros(1, r);
Z(1) = 1;
psialph = poldiv(fliplr(th), fliplr(phialph), r-1);
H = psialph';
[A, Sigma, Xi] = incovma(phi, alpha, th);
