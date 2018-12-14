function [f, e, g, M, AA, Sigma, X] = lkhev(y, Y, Z, T, H, A, Sigma, Xi, d)
%
%    this function computes the residuals and the GLS estimator of a
%    regression model with ARIMA errors using the two-stage Kalman filter.
%    The series can be nonstationary. The state space model is
%
%    x_{t} = T x_{t-1} + H a_{t}
%    y_{t}   = Y \beta + Z x_{t},
%
%    where the initial state vector is
%
%        x_{d+1} = A*\delta + \Xi*c,
%
%    and Var(c) = Sigma. See Gomez and Maravall (1994), "Estimation,
%        Prediction and Interpolation for Nonstationary Series with the
%        Kalman Filter", Journal of the American Statistical Association,
%        89, 611-624. The filter is initialized at time t = d+1, where d is
%        the differencing degree, and the first d observations are stacked
%        to form the \delta vector.
%
%        INPUTS:
%        y: an array containing the input series
%        Y: a matrix containing regression variables
%        Z: the Z matrix
%        T: the T matrix
%        H: the H matrix
%    Sigma: the Sigma matrix
%        A: the A matrix
%       Xi: the Xi matrix
%        d: an integer, the degree of the differencing operator
%
%        OUTPUTS:
%        f: residual vector, whose sum of squares will be minimized
%        e: residual vector for inference
%        g: array containing the regression estimates
%        M: matrix containing the mse of the regression estimates
%       AA: the estimated augmented state vector at the end of filtering
%    Sigma: the Mse of A at the end of filtering
%        X: the design matrix for OLS

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
[r, m] = size(T);
[mY, nY] = size(Y);
n = length(y);
%
% initial covariance matrix
% (note that Xi is a lower triangular matrix)
%
Xis = sparse(Xi);
Sigma = Xis * Sigma * Xis';
%
% initial state vector
%
if d > 0
    if nY > 0
        AA = A * ([Y(1:d, :), y(1:d)]);
    else
        AA = A * y(1:d);
    end
else
    AA = zeros(r, nY+1);
end
%
% Design matrix for OLS
%
X = zeros(n-d, nY+1);
%
% make sparse the last row of T
%
TT = sparse(T(r, :));
%
% Run the Kalman filter
%
[HHp] = mulhkp(H, H);
f = 1;
fc = 0;
K = zeros(r, 1);
KK = K;
for i = d + 1:n
    if nY > 0
        V = [Y(i, :), y(i)] - AA(1, :);
    else
        V = y(i) - AA(1, :);
    end
    F = Sigma(1, 1);
    Fr = sqrt(F);
    f = f * Fr;
    [f, fc] = updatef(f, fc);
    X(i-d, :) = V / Fr;
    KK(1:r-1) = Sigma(2:r, 1);
    KK(r) = TT * Sigma(:, 1);
    K = KK / F;
    wk1 = TT * AA;
    AA(1:r-1, :) = AA(2:r, :);
    AA(r, :) = wk1;
    AA = AA + K * V;
    wk2 = TT * Sigma;
    Sigma(1:r-1, :) = Sigma(2:r, :);
    Sigma(r, :) = wk2;
    wk2 = TT * Sigma';
    Sigma(:, 1:r-1) = Sigma(:, 2:r);
    Sigma(:, r) = wk2';
    [KKp] = mulhkp(K, KK);
    Sigma = Sigma - KKp + HHp;
end
%
% compute the beta estimator, its covariance matrix and the residuals.
% The covariance matrix is not multiplied by sigma^ 2.
%
if nY > 0
    [g, M, e] = bmols(X(:, nY+1), X(:, 1:nY));
else
    e = X;
    g = [];
    M = [];
end
%
% compute determinantal factor
%
% f=f^(1/(n-d));
f = (f^(1 / (n - d))) * (2^(fc / (n - d)));
