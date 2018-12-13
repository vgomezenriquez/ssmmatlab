function [phi, D, DA, ferror] = mecf2mid(Lambda, alpha, betap)
%
% Given the error correction form
%
%       Phi(z) = Lambda(z)*(I-z*I) - Pi*z
%
% corresponding to a polynomial matrix Phi(z) = I + Phi_1*z + Phi_1*z + ...
% + Phi_{p}*z^{p}, where Pi = -Phi(1) = alpha*betap, this function
% obtains the 'differencing' polynomial matrix D(z) = I + D_1*z, were
%
%               D_1 = -betaor*pinv(betaor'*betaor)*betaor',
%
% and betaor is the orthogonal complement of beta = betap', and the
% polynomial matrix phi(z) = I +  phi_1*z + phi_1*z + ...
% + phi_{p-1}*z^{p-1} such that
%
%      Phi(z) =  phi(z)*D(z).
%
% On output, DA = [betaor Idx], where Idx is an index indicating the rows
% of betaor that are linearly independent (=0) and those that are linearly
% dependent (=1).
%
%  Inputst: Lambda : a polynomial matrix of degree p-1
%            alpha : an (s x r) matrix such that Pi = -Phi(1) = alpha*betap
%            betap : an (r x s) matrix
%  Output :    phi : a polynomial matrix of degree p-1
%                D : a polynomial matrix of degree 1
%               DA : an (s x r+1) matrix such that DAf=[betaor Idx]
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

Pi = alpha * betap;
[s, junk] = size(Pi);
beta = betap';
U2(:, :, 1) = eye(s);
U2(:, :, 2) = -beta * pinv(beta'*beta) * beta';
[phi, ferror] = pmatmul(Lambda, U2);
phi(:, :, 2) = phi(:, :, 2) - Pi;
tol = 1d-10;
[phi, gM] = cleanpmat(phi, tol);
%parameterization
[DAb, ferror] = parambeta(beta);
[betaorp, ferror] = m2mor(DAb);
betaor = betaorp';
D(:, :, 1) = eye(s);
D(:, :, 2) = -betaor * pinv(betaor'*betaor) * betaor';
%parameterization
[DA, ferror] = parambeta(betaor);
