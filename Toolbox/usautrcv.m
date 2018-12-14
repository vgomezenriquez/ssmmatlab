function x = usautrcv(xt, m, pfix, pvar, xf)
%**************************************************************************
% This function untransforms the parameters
%
% Auxiliary function called by usa4vcvf.m.
% Reference: ``Estimating Potential Output, Core Inflation
% and the NAIRU as Latent Variables'', by Rafael Domenech
% and Victor Gomez, Journal of Business and Economic Statistics (2006)
%
%    INPUTS:
%       xt : vector containing transformed varying parameters
%        m : number of variables (0,1,2,3 or 4)
%     pfix : array with fixed parameter indices
%     pvar : array with variable parameter indices
%       xf : vector with fixed parameters
%
%    OUTPUT:
%        x : vector with untransformed varying parameters
%
% Copyright (c) 21 July 2003 by Victor Gomez
% Ministerio de Hacienda, Direccion Gral. de Presupuestos,
% Subdireccion Gral. de Analisis y P.E.,
% Alberto Alcocer 2, 1-P, D-34, 28046, Madrid, SPAIN.
% E-mail: VGomez@sepg.minhap.es
%
% The author assumes no responsibility for errors or damage resulting from the use
% of this code. Usage of this code in applications and/or alterations of it should
% be referenced. This code may be redistributed if nothing has been added or
% removed and no money is charged. Positive or negative feedback would be appreciated.
%**************************************************************************

npar = length(pfix) + length(pvar);
xx = zeros(1, npar);
x = zeros(1, length(pvar));
xx(pfix) = xf;
xx(pvar) = xt;

if m == 4
    % four variables
    xx(1:7) = exp(xx(1:7));
    xx(8:9) = (exp(xx(8:9)) - 1) ./ (1 + exp(xx(8:9)));
    r = (exp(xx(10:13)) - 1) ./ (1 + exp(xx(10:13)));
    xx(10:13) = -pafi(r);
    xx(20) = exp(xx(20)) / (1 + exp(xx(20)));
    xx(21) = (0.1963 + 1.0472 * exp(xx(21))) / (1 + exp(xx(21)));
    xx(22:23) = exp(xx(22:23));
    %   xx(22:25)=(1)./(1+exp(xx(22:25)));           % rhos
elseif m == 3
    % three variables
    xx(1:5) = exp(xx(1:5));
    xx(6) = (exp(xx(6)) - 1) ./ (1 + exp(xx(6)));
    r = (exp(xx(7:10)) - 1) ./ (1 + exp(xx(7:10)));
    xx(7:10) = -pafi(r);
    xx(15) = exp(xx(15)) / (1 + exp(xx(15)));
    xx(16) = (0.1963 + 1.0472 * exp(xx(16))) / (1 + exp(xx(16)));
elseif m == 2
    % two variables
    xx(1:3) = exp(xx(1:3));
    r = (exp(xx(4:7)) - 1) ./ (1 + exp(xx(4:7)));
    xx(4:7) = -pafi(r);
    xx(9) = exp(xx(9)) / (1 + exp(xx(9)));
    xx(10) = (0.1963 + 1.0472 * exp(xx(10))) / (1 + exp(xx(10)));
elseif m == 1
    % one variable
    xx(1) = exp(xx(1));
    xx(2) = exp(xx(2)) / (1 + exp(xx(2)));
    xx(3) = (0.1963 + 1.0472 * exp(xx(3))) / (1 + exp(xx(3)));
end
x = xx(pvar);
