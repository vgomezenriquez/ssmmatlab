function str = mhanris2(y, res, x, str)
% PURPOSE: performs the second step of the multivariate Hannan-Rissanen
% method for VARMAX models with restrictions and returns the estimated parameters
%---------------------------------------------------
% USAGE: str = mhanris2(y,res,x,str)
% where:    y      = an (nobs x neqs) matrix of y-vectors
%           res    = an (nobs x neqs) matrix of residuals estimated with a long VARX model
%           x      = matrix of input variables (nobs x nx)
%                   (NOTE: constant vector automatically included)
%           str    = a structure containing the structure of the VARMAX
%---------------------------------------------------
% RETURNS: str = a structure containing the previous structure plus
%                the estimated parameters
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

if nargin ~= 4
    error('wrong # of arguments to mhanris2');
end;

[nobs, neqs] = size(y);
if ~isempty(x)
    [nobs2, nx] = size(x);
    if (nobs2 ~= nobs)
        error('mhanris2: nobs in x-matrix not the same as y-matrix');
    end
else
    nx = 0;
end
[nobs3, nr] = size(res);
if (nobs3 ~= nobs)
    error('mhanris2: nobs in res-matrix not the same as y-matrix');
end;


s = str.s;
m = str.m;
kro = str.kro;

if (s ~= neqs)
    error('mhanris2: s nonequal dim(y_t)');
end
if (m ~= nx)
    error('mhanris2: m nonequal dim(x_t)');
end

str.residv = res;
sigmar = cov(res, 1); %covariance matrix of residuals from the long VAR regression
str.sigmarv = sigmar;


% adjust nobs to feed the lags
nlag = max(kro);
nobse = nobs - nlag;

ylag = glags(y, nlag);
reslag = glags(res, nlag);
vp = y(nlag+1:end, 1:neqs) - res(nlag+1:end, 1:neqs);
if (nx > 0)
    xlag = ltflag(x, nlag);
else
    xlag = [];
end

% form x-matrix, constant is included
% xmat = kron([-ylag -vp reslag xlag ones(nobse,1)],sparse(eye(neqs)));
A = [-ylag, -vp, reslag, xlag, ones(nobse, 1)];
[na, ma] = size(A);
xmat = zeros(na*neqs, ma*neqs);
for i = 1:na
    for j = 1:ma
        aa = repmat(A(i, j), 1, neqs);
        dd = diag(aa);
        xmat((i - 1)*neqs+1:i*neqs, (j - 1)*neqs+1:j*neqs) = dd;
    end
end

% form vector of parameters without restrictions
str = vecparwr(str);
vgam = str.vgam; % the parameters for the constant are at the end of vgam

% form design matrix for the regression by eliminating the columns corresponding
% to the restricted parameters (i.e, those that are zero). Store in bind the indices
% of vgam corresponding to the parameters that are not restricted
zr = [];
bind = [];
for i = 1:length(vgam)
    if isnan(vgam(i))
        zr = [zr, xmat(:, i)];
        bind = [bind; i];
    end
end
str.bind = bind;

[R, p] = chol(sigmar);
L = R';
if (p > 0)
    error('covariance matrix of residuals singular in mhanris2')
end
zrt = []; %[nzr,mzr]=size(zr);
for i = 1:nobse
    zrt = [zrt; L \ zr((i - 1)*neqs+1:i*neqs, :)];
end
ly = L \ y(nlag+1:end, :)';
yt = vec(ly);
[beta, tv] = btval([], [yt, zrt]); % regression
str.beta = beta; % estimated parameters, included the constant at the end
str.tv = tv; % t-values


% insert estimated parameters in vgam
vgams = vgam;
vgamtv = vgam;
[nbind, mbind] = size(bind);
for i = 1:nbind
    vgams(bind(i)) = beta(i);
    vgamtv(bind(i)) = tv(i);
end
% vgams
str.vgams = vgams;
str.vgamtv = vgamtv;

%compute the state space and VARMAX echelon form corresponding to the
% second step
str.noninv2 = 0;
str.nonst2 = 0;


resid2 = zeros(nobse, neqs); % preliminary residuals of the second step regression
resid2v = yt - zrt * beta; % vectorized residuals
for i = 1:nobse
    resid2(i, :) = (L * resid2v((i - 1)*neqs+1:i*neqs))';
end
sigmar2 = cov(resid2, 1);
str.resid2 = resid2;
str.sigmar2 = sigmar2;


if nlag == 0 % the only parameter is the mean
    mu = vgams(end-s+1:end); %parameters for the mean
    mutv = vgamtv(end-s+1:end);
    str.mu = mu;
    str.mutv = mutv;
    return
end

str = param2sse(str);


% test for stationarity and invertibility after the second step.
% if the model is noninvertible, the recursions to obtain the derivatives
% can explode.
% maxmodma=max(abs(eig(str.Fs-str.Ks*str.Hs))); str.maxmodma2=maxmodma;
% maxmodar=max(abs(eig(str.Fs))); str.maxmodar2=maxmodar;
ima = chkstainv(str.Fs-str.Ks*str.Hs); %if ima >1, the model is not invertible
%  report noninvertibility if any
if ima > 1
    str.noninv2 = 1;
end

iar = chkstainv(str.Fs); %if iar >1, the model is not stationary
% report nonstationarity if any
if iar > 1
    str.nonst2 = 1;
end
