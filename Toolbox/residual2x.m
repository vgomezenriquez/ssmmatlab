function [F, e, g, M, A, P, matsis] = residual2x(x, y, Y, s, S, dr, ds, dS, pr, ps, qr, qs, qS)
%
%
%        This function evaluates the residuals for
%        the nonlinear minimization of the sum of
%        squares of an ARIMA model with two possible seasonalities
%
%        INPUTS:
%        x: an array containing model parameters
%        y: an array containing the input series
%        Y: a matrix containing regression variables
%        s:  seasonality
%        S:  second seasonality
%        p:  AR order
%       ps: order of the AR of order s
%        q:  order of the regular MA
%       qs: order of the MA of order s (1 at most)
%       qS: order of the MA of order S (1 at most)
%       dr: order of regular differencing
%       ds: order of differencing of order s
%       dS: order of differencing of order S
%
%        OUTPUTS:
%        F: residual vector, whose sum of squares will be minimized
%        e: residual vector for inference
%        g: array containing the regression estimates
%        M: matrix containing the mse of the regression estimates
%        A: the estimated augmented state vector at the end of filtering
%        P: the Mse of A at the end of filtering
%   matsis: a structure containing the system matrices
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


[phirs, alprsS, thrsS] = arimapol(x, s, S, pr, ps, dr, ds, dS, qr, qs, qS);

d = dr + ds * s + dS * S;

[Z, T, H, A, Sigma, Xi] = arimam(phirs, alprsS, thrsS);

W = [];
G = 0;

matsis.X = Y;
matsis.Z = Z;
matsis.G = G;
matsis.W = W;
matsis.T = T;
matsis.H = H;

% likelihood evaluation, residuals, OLS estimator and initial conditions
%for forecasting
[f, e, g, M, A, P, X] = lkhev(y, Y, Z, T, H, A, Sigma, Xi, d);
F = e * f;
% matsis.OLSMatrix=X;
%OLS residuals
if ~isempty(g)
    nbeta = size(g, 1);
    matsis.olsres = X(:, nbeta+1) - X(:, 1:nbeta) * g;
end

% chb=1.;
% [mz, nalpha]=size(Z);
% Xis=sparse(Xi);
% Sigma=Xis*Sigma*Xis';
% ins=[Sigma A];
% i=[nalpha 0 0 d];
% matsis.ins=ins; matsis.i=i;
% if nY > 0
%  [e,f,g,M,A,P]=scakfle2(y,Y(1:n,:),Z,G,W,T,H,ins,i,chb);
% else
%  [e,f,g,M,A,P]=scakfle2(y,Y,Z,G,W,T,H,ins,i,chb);
% end
% F=e*f;
