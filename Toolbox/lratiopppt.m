function [lagsopt, a, ferror] = lratiopppt(y, x, seas, maxlag, minlag, prt)
% PURPOSE: performs sequential likelihood ratio tests in varmax(p,p,p)
%          models to determine optimal p, starting from p=minlag and
%          proceeding up to p=maxlag. All models are estimated
%          using the same sample size: nobs - maxlag.
%---------------------------------------------------
% USAGE:  [lagsopt,ferror] = lratiopppt(y,x,seas,maxlag,minlag,prt)
% where:    y    = an (nobs x neqs) matrix of y-vectors
%           x    = matrix of input variables (nobs x nx)
%                 (NOTE: constant vector automatically included)
%         seas   = seasonality
%           maxlag = the maximum lag length. If empty on entry, it is
%                    calculated by the program as the order of the VARX
%                    approximation.
%           minlag = the minimum lag length
%           prt = flag for printing
%                    0 = no, 1 = yes
%                    (default = 0)
%---------------------------------------------------
% RETURNS: lagsopt, the optimum number of lags
%                      a = residuals of the VARX approximation (nobs x
%                      neqs)
%---------------------------------------------------
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
%

lagsopt = [];
ferror = 0;

if nargin ~= 6
    disp('wrong # of arguments to lratiopqr');
    ferror = 1;
    return
end;

if maxlag < minlag
    disp('maxlag < minlag in lratiopqr');
    ferror = 2;
    return
end;

[nobs, neqs] = size(y);
if ~isempty(x)
    [mx, nx] = size(x);
    if (mx ~= nobs)
        disp('lratiopqr: nobs in x-matrix not the same as y-matrix');
        ferror = 3;
        return
    end
else
    nx = 0;
end

%First, determination of optimum VARX length
if (nx > 0)
    minp = 5;
else
    minp = 3;
end
if (seas > 1)
    aa = 2.;
else
    aa = 1.25;
end
pt = max(8, seas+minp);
minlags = 0;
maxlags = ceil(max(log(nobs)^aa, pt)); %VARX length is given by this formula
if nx == 0
    [lagsopt, initres] = lratiocr(y, maxlags, minlags, prt);
    res = var_res(y, lagsopt);
else
    [lagsopt, initres] = lratiocrx(y, maxlags, minlags, prt, x);
    res = varx_res(y, lagsopt, x);
end
a = [initres(1:lagsopt, :); res]; %residuals of the VARX approximation
if prt == 1
    fprintf(1, 'estimated order in VARX = %2d\n\n', lagsopt);
end

%the following three lines added on 17-2-2011
if isempty(maxlag)
    maxlag = lagsopt;
end

%identify a varmax(p,p,p)
incr = 1;
lagsopt = lratiocrax(y, maxlag, minlag, incr, prt, a, x);
if prt == 1
    fprintf(1, 'estimated order in VARMAX(p,p,p) = %2d\n\n', lagsopt);
end
