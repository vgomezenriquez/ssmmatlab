function [pr, spr] = predt(n, npr, str, Y, Z, T, H, AA, Sigma, g, M)
%
% this function computes npr forecasts of an ARIMA model and their mse
% using the Kalman filter. The state space model is
%
%    x_{t} = T x_{t-1} + H a_{t}
%    y_{t} = Y_{t} \beta + Z x_{t},
%
%    where he initial state vector is
%
%        x_{d+1} = A*\delta + \Xi*c,
%
% See Gomez and Maravall (1994), "Estimation,
% Prediction and Interpolation for Nonstationary Series with the
% Kalman Filter", Journal of the American Statistical Association,
% 89, 611-624. The filter is initialized at time t = d+1, where d is
% the differencing degree, and the first d observations are stacked
% to form the \delta vector.
%
% Input arguments:
% n  : the series length
% npr: the number of forecasts
% str: estimated standard deviation of the residuals
% Y: matrix containing regression variables
% Z: the Z matrix
% T: the T matrix
% H: the H matrix
% g: array containing the regression estimates
% M: matrix containing the mse of the regression estimates
% AA: the estimated augmented state vector at the end of filtering
% Sigma: the Mse of A at the end of filtering
%
% Output arguments:
% pr : array containing the forecasts
% spr: array containing the mse of the forecasts
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

pr = zeros(npr, 1);
spr = zeros(npr, 1);
[r, m] = size(T);
[mY, nY] = size(Y);
%
% se hace sparse la ultima fila de T
%
TT = sparse(T(r, :));
[HHp] = mulhkp(H, H);
for i = 1:npr
    if nY > 0
        x = AA * [-g; 1];
        pr(i) = x(1) + Y(n+i, :) * g;
        ss = Sigma(1, 1) + (Y(n+i, :) - AA(1, 1:nY)) * M * (Y(n+i, :) - AA(1, 1:nY))';
    else
        x = AA;
        pr(i) = x(1);
        ss = Sigma(1, 1);
    end
    spr(i) = sqrt(ss) * str;
    wk1 = TT * AA;
    AA(1:r-1, :) = AA(2:r, :);
    AA(r, :) = wk1;
    wk2 = TT * Sigma;
    Sigma(1:r-1, :) = Sigma(2:r, :);
    Sigma(r, :) = wk2;
    wk2 = TT * Sigma';
    Sigma(:, 1:r-1) = Sigma(:, 2:r);
    Sigma(:, r) = wk2';
    Sigma = Sigma + HHp;
end
