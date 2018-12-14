
function [A, Sigma, Xi] = incovma(phi, alpha, th)
%
%
%        This function computes the initial covariance matrix, Sigma,
%        and the A matrix for the initial conditions of the Kalman filter.
%        The initial state vector is
%
%        x_{d+1} = A*\delta + \Xi*c,
%
%        where Var(c) = Sigma. See Gomez and Maravall (1994), "Estimation,
%        Prediction and Interpolation for Nonstationary Series with the
%        Kalman Filter", Journal of the American Statistical Association,
%        89, 611-624. The filter is initialized at time t = d+1, where d is
%        the differencing degree, and the first d observations are stacked
%        to form the \delta vector.
%
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
% dimensions for the model
%
p = length(phi) - 1;
q = length(th) - 1;
d = length(alpha) - 1;
r = max(p+d, q+1);
%
% Generate the A_t such that alpha(B)A'_t=0
%
if d > 0
    A = [eye(d); zeros(r, d)];
    for i = d + 1:d + r
        for j = 1:d
            A(i, :) = A(i, :) - alpha(j) * A(i-d+j-1, :);
        end
    end
    A = A(d+1:d+r, :);
else
    A = [];
end
%
% Generate matrix Xi for the initial conditions
%
xi = poldiv(1, fliplr(alpha), r-1);
Xi = eye(r);
for i = 2:r
    Xi(i, i-1:-1:1) = xi(2:i);
end
%
% Initial conditions
%
c = acgf(phi, th, r);
Sigma = zeros(r);
psit = poldiv(fliplr(th), fliplr(phi), r-1);
for i = 1:r
    for j = i:r
        Sigma(i, j) = c(j-i+1) - sum(psit(1:i-1).*psit(j-i+1:j-1));
        Sigma(j, i) = Sigma(i, j);
    end
end
