function [beta, tv, sigmar, covvecbeta, corvecbeta] = multval(y, Y)
%
% This function computes the multiple OLS estimator and the t-values of the
% multivariate linear regression model
%
%  y'_t = Y'_t*beta +  epsilon'_t,       t=1,2,...,n.
%
% Matrices y and Y contain the stack of y'_t and Y'_t, respectively.
%---------------------------------------------------
% USAGE: [beta,tv,sigmar,covvecbeta,corvecbeta]=multval(y,Y)
% where:    y      = an (n x m) matrix of y-vectors
%           Y      = matrix of input variables (n x nY)
%---------------------------------------------------
% RETURNS: beta    = an (nY x m) matrix of regression coefficients
%            tv    = an (nY x m) matrix containing the t-values
%         sigmar   = an (m x m) matrix containing the residual covariance
%                    matrix
%      covvecbeta  = an (nY*m x nY*m) matrix containing the covariance
%                    matrix of vec(beta)
%      corvecbeta  = an (nY*m x nY*m) matrix containing the correlation
%                    matrix of vec(beta)
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
[n, m] = size(y);
[mY, nY] = size(Y);
if n ~= mY
    error('y and Y must have the same number of rows in multval');
end
[beta, M, e] = mulmols(y, Y);
vecbeta = vec(beta);
%the next line changed 9-5-2016
%covariance matrix of residuals
% sigmar=e'*e/(n-nY);
sigmar = e' * e / n;
%end of modification 9-5-2016
%covariance matrix of hat beta
covvecbeta = kron(sigmar, M);
P = diag(sqrt(diag(sigmar)));
rsigmar = (P \ sigmar) / P;
Q = diag(sqrt(diag(M)));
rM = (Q \ M) / Q;
%correlation matrix of hat beta
corvecbeta = kron(rsigmar, rM);
d = sqrt(diag(covvecbeta));
% t-values
tvvecbeta = vecbeta ./ d;
tv = reshape(tvvecbeta, nY, m);
