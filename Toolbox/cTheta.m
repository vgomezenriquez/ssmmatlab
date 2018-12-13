function [Theta, D] = cTheta(n, T, Sigma, phi, th)
%
%    this function computes the lower triangular matrix, Theta, such that
%
%     \Phi*L = \Theta
%
%    using the CKMS recursions corresponing to an ARMA model.
%    The series is assumed to be stationary. The state space model is
%
%    x_{t} = T x_{t-1} + H a_{t}
%    y_{t}   = Z x_{t},
%
%    where Var(x_{1}) = Sigma.
%
%        INPUTS:
%        T: the T matrix
%    Sigma: the Sigma matrix
%      phi: an array containing the AR polynomial
%       th: an array containing the MA polynomial
%
%        OUTPUTS:
%    Theta: lower triangular matrix such that \phi*L = \Theta
%        D: diagonal matrix containing the standar errors
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
Theta = eye(n);
D = zeros(n, 1);
r = size(T, 1);
TT = T(r, :);
idx = find(TT);
TT = TT(idx);
qlength = size(th, 2) - 1;
plength = size(phi, 2) - 1;
m = max(plength, qlength);
%
%extend AR polynomial with zeros if necessary
%
if (m > plength)
    phi = [zeros(1, m-plength), phi];
end
idp = (abs(phi) > 0);
%
% initial covariance matrix is Sigma
% F_1:
F = Sigma(1, 1);
D(1) = F;
K = zeros(r, 1);
L = K;
KK = K;
L(1:r-1) = Sigma(2:r, 1);
L(r) = TT * Sigma(idx, 1);
K = L / F;
%tol: parameter to decide when to pass to steady state. The innovations
%covariance of the filter should converge to one.
tol = 1.d-10;
steady = 0;
dif = 1.d10;
%
% Run the Kalman filter (Morf-Sidhu-Kailath algorithm)
%
for i = 1:n - 1
    ntheta = min(m, n-i);
    theta = zeros(ntheta, 1);
    theta(1) = phi(m) + K(1);
    for j = 2:ntheta
        %   theta(j)=phi(m-j+2:m)*K(1:j-1) + phi(m-j+1) + K(j);
        sum = phi(m-j+1) + K(j);
        for k = m - j + 2:m
            if idp(k)
                sum = sum + phi(k) * K(k-m+j-1);
            end
        end
        theta(j) = sum;
    end
    %  i,theta,pause
    Theta(i+1:i+ntheta, i) = theta;
    if (steady == 0)
        KK(1:r-1) = L(2:r);
        KK(r) = TT * L(idx);
        l1 = L(1);
        l1f = l1 / F;
        L = KK - l1 * K;
        F = F - l1f * l1;
        K = K - (l1f) * L / F;
        if (dif <= tol)
            steady = 1; %no more updates, steady state instead
        end
        dif = abs(F-1.d0);
    end
    D(i+1) = F;
end
