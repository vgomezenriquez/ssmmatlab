function x0 = cinest(y, Y, parm, est, ols, a, prt, fid)
%
% function to estimate initial parameter values in an ARIMA model
%
%
%     INPUTS:
%     y       : input series
%     Y       : matrix of regression variables
%     parm    : a structure, containing the ARIMA specification. It should
%               have al least the following fields:
%          .s : seasonality
%          .S : second seasonality
%          .p : degree of regular AR polynomial
%          .d : degree of regular differencing
%          .q : degree of regular MA polynomial
%          .ps: degree of seasonal AR polynomial
%          .ds: degree of seasonal differencing
%          .qs: degree of seasonal MA polynomial
%          .dS: degree of second seasonal differencing
%          .qS: degree of second seasonal MA polynomial
%     est = 1 : estimation of regression coefficients
%        = 0  : no estimation of regression coefficients
%        ols  : = 1, perform OLS, = 0, use the Durbin Levinson algorithm in the HR
%               method
%     a       : an integer, the degree of the AR approximation in the first
%               step of the Hanna-Rissanen method.
%     prt     : = 1, print results in an external file, = 0, do not print
%     fid     : an integer, corresponding to an external file
%
%  OUTPUTS:
%      x0     : an array containg the initial estimated values
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

if ~isstruct(parm)
    error('cinest: requires a parameter structure');
end;
s = parm.s;
S = parm.S;
dr = parm.dr;
ds = parm.ds;
dS = parm.dS;
p = parm.p;
ps = parm.ps;
q = parm.q;
qs = parm.qs;
qS = parm.qS;
[yd, beta] = diffest(y, Y, s, S, dr, ds, dS, est);
x0 = inest(yd, beta, s, S, p, ps, q, qs, qS, ols, a);
chk = chkroots(x0, p, ps, q, qs, qS); %check stationarity and invertibility
if chk == 1, x0 = invroots(x0, p, ps, q, qs, qS);
end
if prt == 1
    fprintf(fid, 'Initial parameter values:\n');
    fprintf(fid, '%12.4f', x0);
    fprintf(fid, '\n\n');
end
