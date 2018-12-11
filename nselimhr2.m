function [mp, minc, mint] = nselimhr2(y, x, str)
% PURPOSE: eliminates nonsignificant parameters in the second stage of the
% Hannan-Rissanen method for VARMAX models with restrictions
%---------------------------------------------------
% USAGE:  [mp, minc, mint] = nselimhr2(y,x,str)
% where:    y      = an (nobs x neqs) matrix of y-vectors
%                x      = matrix of input variables (nobs x nx)
%                   (NOTE: constant vector automatically included)
%               str    = a structure containing the structure of the VARMAX   model
%---------------------------------------------------
% RETURNS:
%           minc  = a 1 x 3 array containing the index in the array of the eliminated
%                         parameter
%             minc = 'ph', 'th' or 'ga', referring to AR, MA or X part to
%                         which the parameter belongs
%             mint =  the t-value of the eliminated parameter
%---------------------------------------------------
% Copyright (c) 21 July 2003 by Victor Gomez
% Ministerio de Hacienda, Direccion Gral. de Presupuestos,
% Subdireccion Gral. de Analisis y P.E.,
% Alberto Alcocer 2, 1-P, D-34, 28046, Madrid, SPAIN.
% Phone : +34-915835439
% E-mail: vgomez@sgpg.meh.es
%
% The author assumes no responsibility for errors or damage resulting from the use
% of this code. Usage of this code in applications and/or alterations of it should
% be referenced. This code may be redistributed if nothing has been added or
% removed and no money is charged. Positive or negative feedback would be appreciated.
%

if nargin ~= 3
    error('wrong # of arguments nselimhr');
end;

[nobs, neqs] = size(y);
if ~isempty(x)
    [nobs2, nx] = size(x);
    if (nobs2 ~= nobs)
        error('nselimhr2: nobs in x-matrix not the same as y-matrix');
    end
else
    nx = 0;
end

s = str.s;
m = str.m;
phi = str.phi;
theta = str.theta;
gamma = str.gamma;
kro = str.kro;
phitv = str.phitv;
thetatv = str.thetatv;
gammatv = str.gammatv;

if (s ~= neqs)
    error('nselimhr2: s nonequal dim(y_t)');
end
if (m ~= nx)
    error('nselimhr2: m nonequal dim(x_t)');
end

minc = [0, 0, 0];
mint = 1.d15;
mp = [];
[np, lp, rp] = size(phi);
[nq, lq, rq] = size(theta);
[ng, lg, rg] = size(gamma);

it = np;
jt = lp;
mt = rp;
for i = 1:it
    for j = 1:jt
        for k = 2:mt
            tval = abs(phitv(i, j, k));
            if (isnan(phi(i, j, k)) && (tval < mint))
                nkro = 0;
                if k == kro(i) + 1 %enforce row degrees
                    nkro = 1;
                    for jj = 1:jt
                        if (jj ~= j) && (isnan(phi(i, jj, k)))
                            nkro = 0;
                        end
                        if (nx > 0)
                            if (jj <= lq)
                                if (isnan(theta(i, jj, k)))
                                    nkro = 0;
                                end
                            end
                            if (jj <= lg)
                                if (isnan(gamma(i, jj, k)))
                                    nkro = 0;
                                end
                            end
                        else
                            if (jj <= lq)
                                if (isnan(theta(i, jj, k)))
                                    nkro = 0;
                                end
                            end
                        end
                    end
                end
                if nkro == 0
                    minc = [i, j, k];
                    mint = tval;
                    mp = 'ph';
                end
            end
        end
    end
end

it = nq;
jt = lq;
mt = rq;
for i = 1:it
    for j = 1:jt
        for k = 1:mt
            tval = abs(thetatv(i, j, k));
            if (isnan(theta(i, j, k)) && (tval < mint))
                nkro = 0;
                if k == kro(i) + 1 %enforce row degrees
                    nkro = 1;
                    for jj = 1:jt
                        if (jj ~= j) && (isnan(theta(i, jj, k)))
                            nkro = 0;
                        end
                        if (nx > 0)
                            if (jj <= lp)
                                if (isnan(phi(i, jj, k)))
                                    nkro = 0;
                                end
                            end
                            if (jj <= lg)
                                if (isnan(gamma(i, jj, k)))
                                    nkro = 0;
                                end
                            end
                        else
                            if (jj <= lp)
                                if (isnan(phi(i, jj, k)))
                                    nkro = 0;
                                end
                            end
                        end
                    end
                end
                if nkro == 0
                    minc = [i, j, k];
                    mint = tval;
                    mp = 'th';
                end
            end
        end
    end
end

it = ng;
jt = lg;
mt = rg;
for i = 1:it
    for j = 1:jt
        for k = 1:mt
            tval = abs(gammatv(i, j, k));
            if (isnan(gamma(i, j, k)) && (tval < mint))
                nkro = 0;
                if k == kro(i) + 1 %enforce row degrees
                    nkro = 1;
                    for jj = 1:jt
                        if (jj ~= j) && (isnan(gamma(i, jj, k)))
                            nkro = 0;
                        end
                        if (jj <= lp)
                            if (isnan(phi(i, jj, k)))
                                nkro = 0;
                            end
                        end
                        if (jj <= lq)
                            if (isnan(theta(i, jj, k)))
                                nkro = 0;
                            end
                        end
                    end
                end
                if nkro == 0
                    minc = [i, j, k];
                    mint = tval;
                    mp = 'ga';
                end
            end
        end
    end
end
