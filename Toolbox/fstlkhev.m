function [f, e, g, M] = fstlkhev(y, Y, T, Sigma, chb, nmiss)
%
%    this function computes the residuals using the CKMS recursions
%    corresponing to an ARIMA model, possibly with regression variables.
%    The series is assumed to be stationary. The state space model is
%
%    x_{t} = T x_{t-1} + H a_{t}
%    y_{t}   = Y \beta + Z x_{t},
%
%    where Var(x_{1}) = Sigma.
%
%        INPUTS:
%        y: an array containing the input series
%        Y: a matrix containing regression variables
%        Z: the Z matrix
%        T: the T matrix
%        H: the H matrix
%    Sigma: the Sigma matrix
%      chb: = 1, compute regression estimates, = 0, do not compute reg.
%           estimates.
%
%        OUTPUTS:
%        f: residual vector, whose sum of squares will be minimized
%        e: residual vector for inference
%        g: array containing the regression estimates
%        M: matrix containing the mse of the regression estimates
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
if nargin < 6
    nmiss = 0;
end
r = size(T, 1);
nY = size(Y, 2);
n = length(y);
%
% make sparse the last row of T
%
% TT=sparse(T(r,:));
TT = T(r, :);
idx = find(TT);
TT = TT(idx);
%
% initial covariance matrix is Sigma
% F_1:
F = Sigma(1, 1);
% L_1, K_1:
K = zeros(r, 1);
L = K;
KK = K;
L(1:r-1) = Sigma(2:r, 1);
% L(r)=TT*Sigma(:,1);
L(r) = TT * Sigma(idx, 1);
K = L / F;
%
% initial state vector
%
AA = zeros(r, nY+1);
%
% Design matrix for OLS
%
X = zeros(n, nY+1);
%tol: parameter to decide when to pass to steady state. The innovations
%covariance of the filter should converge to one.
tol = 1.d-5;
steady = 0;
dif = 1.d10;
%
% Run the Kalman filter (Morf-Sidhu-Kailath algorithm)
%
f = 1;
fc = 0;
for i = 1:n
    if nY > 0
        V = [Y(i, :), y(i)] - AA(1, :);
    else
        V = y(i) - AA(1, :);
    end
    %     wk1=TT*AA;
    wk1 = TT * AA(idx, :);
    AA(1:r-1, :) = AA(2:r, :);
    AA(r, :) = wk1;
    AA = AA + K * V;
    if (steady == 0)
        Fr = sqrt(F);
        f = f * Fr;
        [f, fc] = updatef(f, fc);
        %      X(i,:)=V/Fr;
        KK(1:r-1) = L(2:r);
        %      KK(r)=TT*L;
        KK(r) = TT * L(idx);
        l1 = L(1);
        l1f = l1 / F;
        L = KK - l1 * K;
        F = F - l1f * l1;
        K = K - (l1f) * L / F;
        if (dif <= tol)
            steady = 1; %no more updates, steady state instead
        end
        dif = abs(Fr-1.d0);
    end
    X(i, :) = V / Fr;
end
%
% compute the residuals
%
if nY > 0
    [Q, R] = qr(X(:, 1:nY));
    qy = Q' * X(:, nY+1);
    e = qy(nY+1:n);
    if chb == 1
        g = R(1:nY, :) \ qy(1:nY);
        M = pinv(R(1:nY, :));
        M = M * M';
    else
        g = [];
        M = [];
    end
else
    e = X;
    g = [];
    M = [];
end
%
% compute determinantal factor
%
%f=f^(1/n);
if (nmiss > 0)
    %determinantal correction for missing values
    LL = R(1:nmiss, 1:nmiss);
    f = f * abs(prod(diag(LL)));
    [f, fc] = updatef(f, fc);
end
f = (f^(1 / n)) * (2^(fc / n));
