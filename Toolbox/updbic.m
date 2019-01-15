function [bicm, oparm] = updbic(yd, beta, s, S, p, ps, q, qs, qS, ols, a, bicm, oparm)
%
% this function updates the BIC criterion of an ARMA model after checking
% for stationarity and invertibility
%
% Input arguments:
% yd   : an (n x m) matrix containing the series, yd(:,1), and an (n x m-1)
%        matrix of regression variables if m > 1.
% beta : an m-1 vector containing the OLS estimators if m > 1, empty if m =
%        1
% s    : seasonality
% S    : second seasonality
% p    : degree of AR polynomial
% ps   : degree of AR seasonal polynomial
% q    : degree of MA polynomial
% qs   : degree of MA seasonal polynomial
% qS   : degree of MA second seasonal polynomial
% ols  : = 1, perform OLS, = 0, use the Durbin Levinson algorithm
% a    : an integer, the degree of the AR approximation in the first step
%        of the Hanna-Rissanen method.
% bicm : the previous bic criterion
% oparm: a structure where
% .s:  seasonality
% .S:  second seasonality
% .qS: order of the MA of order S (1 at most)
% .dr: order of regular differencing
% .ds: order of differencing of order s
% .dS: order of differencing of order S
% .p:  AR order
% .ps: order of the AR of order s
% .q:  order of the regular MA
% .qs: order of the MA of order s (1 at most)
% .pvar:  array containing the indices of variable parameters
% .pfix:  array containing the indices of fixed parameters
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

[nd, junk] = size(yd);
chk = 0;
if p + ps + q + qs > 0
    x = inest(yd, beta, s, S, p, ps, q, qs, qS, ols, a);
    chk = chkroots(x, p, ps, q, qs, qS); %check stationarity and invertibility
else
    x = [];
end
if chk == 0
    bic = cbic(x, yd, nd, s, p, ps, q, qs);
    %  rbic=abs(1-bic/bicm);
    %  if (bic < bicm)
    %   if  (((rbic > .001) & ((rbic > 0.003) | (ps+qs <= oparm.ps+oparm.qs))) |...
    %         ((rbic <= .001) & (ps+qs <= oparm.ps+oparm.qs)))
    %   if  ((rbic > 0.003) | ((rbic <= .003) & (ps+qs <= oparm.ps+oparm.qs)))
    bicm = bic;
    oparm.p = p;
    oparm.q = q;
    oparm.ps = ps;
    oparm.qs = qs;
    %   end
    %  end
end
