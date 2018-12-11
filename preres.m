function [Z, T, H, A, Sigma] = preres(x, s, pr, ps, qr, qs)
%
%
%  This function computes the system matrices corresponding to a stationary
%  ARMA model. The state space model is
%
%    x_{t} = T x_{t-1} + H a_{t}
%    y_{t} = Z x_{t},
%
%    where Var(x_{1}) = Sigma. Given an ARIMA model,the initial state
%    vector is
%
%        x_{d+1} = A*\delta + \Xi*c,
%
% and Var(c) = Sigma. See Gomez and Maravall (1994), "Estimation,
% Prediction and Interpolation for Nonstationary Series with the
% Kalman Filter", Journal of the American Statistical Association,
% 89, 611-624. The filter is initialized at time t = d+1, where d is
% the differencing degree, and the first d observations are stacked
% to form the \delta vector. In this case, because d = 0, x_1 = c;

%
%        INPUTS:
%        x: array containing model parameters
%        s:  number of seasons
%       pr:  AR order
%       ps: order of the AR of order s
%       qr:  order of the regular MA
%       qs: order of the MA of order s
%
%        OUTPUTS:
%        Z: the Z matrix
%        T: the T matrix
%        H: the H matrix
%        A: the empty matrix because the series is stationary
%    Sigma: the Sigma matrix
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

alprs = 1;
phi = 1;
phis = 1;
th = 1;
ths = 1;
if pr > 0
    phi = conv(phi, [fliplr(x(1:pr)), 1]);
end
prs = pr + ps;
if qr > 0
    th = conv(th, [fliplr(x(prs+1:prs+qr)), 1]);
end
if s > 1
    if ps > 0
        phis = zeros(1, s+1);
        phis(1) = x(pr+ps);
        phis(s+1) = 1;
    end
    if qs > 0
        ths = zeros(1, s+1);
        ths(1) = x(prs+qr+qs);
        ths(s+1) = 1;
    end
end
phirs = conv(phi, phis);
thrs = conv(th, ths);


[Z, T, H, A, Sigma] = arimam(phirs, alprs, thrs);
