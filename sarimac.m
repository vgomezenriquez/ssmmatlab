function x = sarimac(p, ps, q, qs, phir, phis, thr, ths)
%**************************************************************************
% Auxiliary function called in arimasimul_d.m to set the ARIMA coefficients
%
%  INPUTS:
%           p:  AR order
%          ps:  order of the AR of order s
%           q:  order of the regular MA
%          qs:  order of the MA of order s (1 at most)
%     phir   : an array containing the regular AR polynomial
%     phis   : an array containing the seasonal AR polynomial
%     thr    : an array containing the regular MA polynomial
%     thr    : an array containing the seasonal MA polynomial
%
%  OUTPUTS:
%           x:  an array containing the ARIMA parameter values
%
% Copyright (c) 21 July 2003 by Victor Gomez
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
%*************************************************************************

x = zeros(1, p+ps+q+qs);
if p > 0
    pp = length(phir) - 1;
    if pp ~= p, error('phir has a length different from p');
    end
    x(1:p) = fliplr(phir(1:p));
end
if ps > 0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % seasonal AR polynomial
    x(p+1:p+ps) = fliplr(phis(1:ps));
end
if q > 0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % regular MA polynomial
    %
    qq = length(thr) - 1;
    if qq ~= q, error('thr has a length different from q');
    end
    x(p+ps+1:p+ps+q) = fliplr(thr(1:q));
end
if qs > 0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % seasonal MA polynomial
    %
    x(p+ps+q+1:p+ps+q+qs) = fliplr(ths(1:qs));
end
