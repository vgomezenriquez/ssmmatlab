function [Pi, Lambda, alpha, betap, ferror] = mid2mecf(phi, D, DAf)
%
% Given a polynomial matrix of the form phi(z)*D(z), where phi(z) = I +
% phi_1*z + phi_1*z + ... + phi_{p}*z^{p} is an autoregressive polynomial
% matrix and D(z) = I + D_1*z is a 'differencing' polynomial matrix, this
% function obtains the error correction form such that
%
%       phi(z)*D(z) = Lambda(z)*(I-z*I) - Pi*z,
%
% where Lambda(z) = I + Lambda_1*z + ... + Lambda_{p-1}*z^{p-1} and Pi =
% -phi(1)*D(1). It is assumed that
%
%               D_1 = -betaor*pinv(betaor'*betaor)*betaor',
%
% where DAf = [betaor Idx] and Idx is an index indicating the rows of
% betaor that or linearly independent.
%
% Inputs  :    phi: a polynomial matrix of degree p
%                D: a polynomial matrix of degree 1
%              DAf: an (s x r+1) matrix such that DAf=[betaor Idx]
%  Output :    Pi : a matrix such that Pi = -phi(1)*D(1) = alpha*betap
%          Lambda : a polynomial matrix of degree p-1
%           alpha : an (s x r) matrix
%           betap : an (r x s) matrix
%
% Copyright (c) 21 July 2014 by Victor Gomez
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

[ny, my, ky] = size(phi);
phibar = pmatmul(phi, D);
[npb, mpb, kpb] = size(phibar);
Pi = -sum(phi, 3) * sum(D, 3);
Aux = phibar;
Aux(:, :, 2) = Aux(:, :, 2) + Pi;
Dif(:, :, 1) = eye(my);
Dif(:, :, 2) = -eye(my);
[Lambda, ierror] = prtransfer(Dif, Aux, kpb-1);
tol = 1d-10;
[Lambda, gM] = cleanpmat(Lambda, tol);
[betap, ferror] = m2mor(DAf);
alpha = Pi / betap;
