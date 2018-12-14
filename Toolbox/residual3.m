function [F, e] = residual3(x, y, Y, s, pr, ps, qr, qs)
%
%
%        This function evaluates the residuals for
%        the nonlinear minimization of the sum of
%        squares of an ARMA model using the CKMS recursions
%
% Input arguments:
% x: array containing model parameters
% y: vector containing the data
% Y: matrix containing regression variables
% s:  number of seasons
% pr:  AR order
% ps: order of the AR of order s
% qr:  order of the regular MA
% qs: order of the MA of order s
%
% Output arguments:
%  F: residual vector, whose sum of squares will be minimized
%  e: residual vector for inference
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

[~, T, ~, ~, Sigma] = preres(x, s, pr, ps, qr, qs);
%perform fast likelihood evaluation
[F, e] = fstlkhev(y, Y, T, Sigma);
