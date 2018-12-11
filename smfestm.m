function [F, xv] = smfestm(xv, y, pfix, pvar, xf, chb, models)
%*************************************************************************
% This function transforms model parameters so that their values lie in the
% feasable region and evaluates the residuals for the nonlinear
% minimization of the sum of squares of a canonical structural model
%
%   INPUTS:
%      xv : vector with parameters to be estimated
%       y : data vector
%       s : frequency of the data
%    pfix : array with fixed parameter indices
%    pvar : array with variable parameter indices
%      xf : vector with fixed parameters
%     chb = 1 : compute beta estimator hb and MSE of hb
%         = 0 : do not compute hb and MSE of hb
%  models : structure containing model information
%
%  OUTPUTS:
%       F : residual vector multiplied with factor f given by scakfle2;
%           used in minimization of nonlinear sum of squares
%      xv : vector with parameters to be estimated that has been tranformed
%           if necessary
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
%**************************************************************************

arp = models.arp;
cycle = models.cycle;
seas = models.seas;
if ~isempty(seas)
    N = length(seas);
else
    N = 0;
end
if (arp + cycle) > 0
    npar = length(pfix) + length(pvar);
    x = zeros(1, npar);
    x(pfix) = xf;
    x(pvar) = xv;
    if arp > 0
        xar = x(end-arp+1:end);
        chk = chkroots(xar, arp, 0, 0, 0, 0);
        if chk == 1
            x(end-arp+1:end) = invroots(xar, arp, 0, 0, 0, 0);
            xv = x(pvar);
        end
    end
    if cycle > 0
        stord = models.stord;
        xp = zeros(1, max(stord));
        xp(stord) = x;
        xrho = xp(6+N);
        if (xrho < 0.) | (xrho > 1.)
            xp(6+N) = exp(xrho) / (1. + exp(xrho));
            x = xp(stord);
            xv = x(pvar);
        end
        xlambda = xp(7+N);
        xl1 = models.xl1;
        xl2 = models.xl2;
        if (xlambda < xl1) | (xlambda > xl2)
            xp(7+N) = (xl1 + xl2 * exp(xlambda)) / (1. + exp(xlambda));
            x = xp(stord);
            xv = x(pvar);
        end
    end
end

[F, e, g, M, Pevf] = smfunm(xv, y, pfix, pvar, xf, chb, models);
