function [resid, sigmar] = compresde0(y, x, str)
% PURPOSE: given a structure, it computes the
% model residuals and their covariance matrix using the difference
% equation, starting with zeros.
%---------------------------------------------------
% USAGE: [resid2,sigmar2]=compresde0(y,x,str)
% where:    y      = an (nobs x neqs) matrix of y-vectors
%           x      = matrix of input variables (nobs x nx)
%           str    = a structure containing the model information
%---------------------------------------------------
% RETURNS: resid  = the residuals
%          sigmar = the residuals covariance matrix
%---------------------------------------------------
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

[nobs, neqs] = size(y);
[nobs2, nx] = size(x);
nlag = max(str.kro);

mu = str.mu; %mu is the constant, mu=phi(1)*med, where med is the mean
if isfield(str, 'phis')
    phi1 = sum(str.phis, 3); %this is phis(1), where phis is the echelon form AR poly.
else
    phi1 = eye(neqs);
end
if ~any(any(phi1-eye(neqs)))
    med = mu;
else
    med = phi1 \ mu;
end
yc = y - repmat(med', nobs, 1); %center variables


if nlag == 0
    resid = yc;
    sigmar = cov(resid, 1);
    return
end


%augment resid, x and y with zeros
npl = nobs + nlag;
nlp1 = nlag + 1;
resida = zeros(npl, neqs);
ya = zeros(npl, neqs);
ya(nlp1:npl, :) = yc;
if (nx > 0)
    xa = zeros(npl, nx);
    xa(nlp1:npl, :) = x;
end

thetast = str.thetast;
phist = str.phist;
gammast = str.gammast;
Th = [];
Ph = [];
Gm = [];
for i = 2:nlp1
    Th = [Th, thetast(:, :, i)];
    Ph = [Ph, phist(:, :, i)];
    if (nx > 0)
        Gm = [Gm, gammast(:, :, i)];
    end
end

for i = nlp1:npl
    Rl = [];
    Yl = [];
    Xl = [];
    for j = i - 1:-1:i - nlag
        Rl = [Rl, resida(j, :)];
        Yl = [Yl, ya(j, :)];
        if (nx > 0)
            Xl = [Xl, xa(j, :)];
        end
    end
    suma = ya(i, :) + Yl * Ph' - Rl * Th'; % constant is included
    if (nx > 0)
        suma = suma - xa(i, :) * gammast(:, :, 1)' - Xl * Gm';
    end
    resida(i, :) = suma;
end
resid = resida(nlp1:npl, :);
%covariance matrix of residuals from the second step regression
sigmar = cov(resid, 1);
