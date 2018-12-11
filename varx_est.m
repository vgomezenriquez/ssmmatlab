function res = varx_est(y, nlag, x, test, xx)
% PURPOSE: performs vector autogressive with exogenous inputs (VARX)
%          estimation and returns a structure
%---------------------------------------------------
% USAGE: res = varx_est(y,nlag,x,test,xx)
% where:    y    = an (nobs x neqs) matrix of y-vectors
%           nlag = the lag length
%           x    = matrix of input variables (nobs x nx)
%                 (NOTE: constant vector automatically included)
%           test = a logical variable to perform additional tests
%           xx   = optional matrix of variables (nobs x nxx)
%---------------------------------------------------
% RETURNS: a structure containing the following fields
%          .resid   = residuals
%          .phi       = VARX matrix polynomials corresponding to the outputs.
%                          Signs are those of Box-Jenkins  (I-phi_1*z-phi_2*z^2 -
%                           ...-phi_p*z^p)
%          .phitv    = matrix polynomials containing the t-values and
%                          corresponding to the output VARX matrix polynomials
%          .phix      = VARX matrix polynomials corresponding to the inputs.
%                           Signs are those of Box-Jenkins  (phix_0 -phix_1*z-
%                           -phix_2*z^2 -...-phix_p*z^p)
%          .phixtv   = matrix polynomials containing the t-values and
%                          corresponding to the input VARX matrix polynomials
%          .const   = vector containing the estimated constant
%          .consttv = vector containing the t-values of the estimated
%                            constant
%          .betava  = matrix containing the estimated regression
%                        coefficients, beta
%          .tvvar     = t-values of beta
%          .sigmar   = covariance matrix of residuals
%          .covvecbeta = covariance matrix of vec(beta)
%          .corvecbeta = correlation matrix of vec(beta)
%          .dusigmar   = determinant of the maximum likelihood estimator of
%                        the covariance matrix of residuals
%          .llkhd      = log-likelihood
%          .aic        = aic
%          .bic        = bic
%          .ssr(j)     = Sum-of-squares residuals for each equation. Only
%                        if test = 1
%          .Rsqr(j)    = R^2 for each equation. Only if test = 1
%          .SEeq(j)    = Standard error of regression  for each equation.
%                        Only if test = 1
%          .Fstat(j)   = F-statistic for each equation. Only if test = 1
%          .llkhdeq(j) = log-likelihood for each equation. Only if test = 1
%          .aiceq(j)   = aic for each equation. Only if test = 1
%          .biceq(j)   = bic for each equation. Only if test = 1
%          .ftest(i,j) = Granger F-tests, only if test = 1
%          .fprob(i,j) = Granger marginal probabilities, only if test = 1
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

if isempty(x)
    if nargin == 5
        res = var_est(y, nlag, test, xx);
    elseif nargin == 4
        res = var_est(y, nlag, test);
    elseif nargin == 3
        res = var_est(y, nlag);
    end
    return
end

[nobs, neqs] = size(y);

mtest = 0;
if nargin >= 4
    mtest = test;
end

nxx = 0;
if nargin == 5
    [nobs2xx, nxx] = size(xx);
    if (nobs2xx ~= nobs)
        error('varx_est: nobs in xx-matrix not the same as y-matrix');
    end;
end;

[nobs2, nx] = size(x);
if (nobs2 ~= nobs)
    error('varx_est: nobs in x-matrix not the same as y-matrix');
end;

% adjust nobs to feed the lags
nobse = nobs - nlag;
dn = double(nobse);

ylag = glags(y, nlag);
xlag = ltflag(x, nlag);

% form x-matrix
if nxx > 0
    xmat = [ylag, xlag, xx(nlag+1:nobs, :), ones(nobse, 1)];
else
    xmat = [ylag, xlag, ones(nobse, 1)];
end;


[mxm, nxm] = size(xmat);

res.resid = zeros(nobse, neqs);
res.betavar = zeros(nxm, neqs);
res.tvvar = zeros(nxm, neqs);

[beta, tv, sigmar, covvecbeta, corvecbeta] = multval(y(nlag+1:nobs, :), xmat);
res.resid = y(nlag+1:nobs, :) - xmat * beta;
res.betavar = beta;
res.tvvar = tv;
res.sigmar = sigmar;
res.covvecbeta = covvecbeta;
res.corvecbeta = corvecbeta;

%VARX coefficient matrices. Signs are those of Box-Jenkins (I-phi_1*z
%-phi_2*z^2 - ...-phi_p*z^p)
phi(:, :, nlag+1) = zeros(neqs);
phi(:, :, 1) = eye(neqs);
phitv = phi;
phix(:, :, nlag+1) = zeros(neqs, nx);
phixtv = phix;
for i = 1:nlag
    ii = (i - 1) * neqs + 1:i * neqs;
    phi(:, :, i+1) = beta(ii, :)';
    phitv(:, :, i+1) = tv(ii, :)';
end
orx = nlag * neqs;
for i = 1:nlag + 1
    ii = orx + (i - 1) * nx + 1:orx + i * nx;
    phix(:, :, i) = beta(ii, :)';
    phixtv(:, :, i) = tv(ii, :)';
end
const = beta(end, :)';
consttv = tv(end, :)';
res.phi = phi;
res.phix = phix;
res.const = const;
res.phitv = phitv;
res.phixtv = phixtv;
res.consttv = consttv;


%loglikelihood
usigmar = sigmar * (dn - double(nxm)) / dn; %maximum likelihood estimator of
%residual covariance matrix
dusigmar = det(usigmar);
dneqs = double(neqs);
l = -(dn / 2) * (dneqs * (1 + log(2*pi)) + log(dusigmar));
npar = dneqs * double(nxm);
%aic
aic = -2 * l / dn + 2 * npar / dn;
%bic
bic = -2 * l / dn + log(dn) * npar / dn;

res.dusigmar = dusigmar;
res.llkhd = l;
res.aic = aic;
res.bic = bic;

if mtest == 1
    
    % additional tests
    res.ssr = zeros(1, neqs);
    res.Rsqr = zeros(1, neqs);
    res.SEeq = zeros(1, neqs);
    res.Fstat = zeros(1, neqs);
    res.llkhdeq = zeros(1, neqs);
    res.aiceq = zeros(1, neqs);
    res.biceq = zeros(1, neqs);
    res.ftest = zeros(neqs, neqs);
    res.fprob = zeros(neqs, neqs);
    nparm = neqs * nlag + 1 + nx;
    dnparm = double(nparm);
    df = dn - dnparm;
    
    % pull out each y-vector and obtain statistics for each equation
    for j = 1:neqs;
        yvec = y(nlag+1:nobs, j);
        resid = res.resid(:, j);
        sigu = resid' * resid; %Sum-of-squared residuals
        res.ssr(:, j) = sigu;
        yvecm = yvec - mean(yvec);
        rsqr1 = sigu;
        rsqr2 = yvecm' * yvecm;
        Rsqr = 1.0 - rsqr1 / rsqr2; %R-squared
        res.Rsqr(:, j) = Rsqr;
        SEeq = sqrt(sigu/df); %Standard Error of regression
        res.SEeq(:, j) = SEeq;
        Fstat = (Rsqr / double(nparm-1)) / ((1 - Rsqr) / df); %F-statistic
        res.Fstat(:, j) = Fstat;
        l = -(dn / 2) * (1 + log(2*pi) + log(sigu/dn));
        res.llkhdeq(:, j) = l;
        aic = -2 * l / dn + 2 * dnparm / dn;
        res.aiceq(:, j) = aic;
        bic = -2 * l / dn + log(dn) * dnparm / dn;
        res.biceq(:, j) = bic;
        % form x matrices for joint F-tests
        % exclude each variable from the model sequentially
        for r = 1:neqs;
            xtmp = [];
            for s = 1:neqs;
                if s ~= r
                    xlag = glags(y(:, s), nlag);
                    xtmp = [xtmp, xlag];
                end;
            end;
            % we have an xtmp matrix that excludes 1 variable
            % add deterministic variables (if any) and constant term
            if nx > 0
                xtmp = [xtmp, x(nlag+1:nobs, :), ones(nobse, 1)];
            else
                xtmp = [xtmp, ones(nobse, 1)];
            end;
            % get ols residual vector
            etmp = rbols(yvec, xtmp);
            sigr = etmp' * etmp;
            % joint F-test for variable r
            ftest(r, 1) = ((sigr - sigu) / nlag) / (sigu / df);
        end;
        res.ftest(:, j) = ftest;
        pf = fdis_cdf(ftest, nlag, df);
        res.fprob(:, j) = ones(neqs, 1) - pf;
    end;
    % end of loop over equations
end
