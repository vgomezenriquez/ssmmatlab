function [fjac, g] = conmedfjac(beta, y, x, str)
% PURPOSE: this function computes the jacobian and the gradient
% corresponding to a VARMAX model
%---------------------------------------------------
% USAGE: [fjac,g]=conmedfjac(beta,y,x,str)
% where:    beta   = a (1 x nparma) vector of parameters
%           y      = an (nobs x neqs) matrix of y-vectors
%           x      = matrix of input variables (nobs x nx)
%           str    = a structure containing the model information
%---------------------------------------------------
% RETURNS: fjac  = the jacobian
%          g     = the gradient
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

kro = str.kro;
nlag = max(kro);
vgams = str.vgam;
bind = str.bind;
[nbind, mbind] = size(bind);
for i = 1:nbind
    vgams(bind(i)) = beta(i);
end
str.vgams = vgams;
str = param2sse(str);

%check stationarity and invertibility. If necessary, change parameters.
% iar=chkstainv(str.Fs);    %if iar >1, the model is not stationary
% if iar > 1
% %  fprintf(1,'model nonstationary, iar = %2d\n',iar);
%  vgam=enfstab(str,'phi  ');
%  str.vgams=vgam; str=param2sse(str);
% end
ima = chkstainv(str.Fs-str.Ks*str.Hs); %if ima >1, the model is not invertible
if ima > 1
    %  fprintf(1,'model noninvertible, ima = %2d\n',ima);
    vgam = enfstab(str, 'theta');
    str.vgams = vgam;
    str = param2sse(str);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% new second step residuals. We compute them starting with zeros. Do not
% use the Kalman filter at this stage.
[resid2, sigmar2] = compresde0(y, x, str); % constant is passed in str.mu
% plot(resid2)
% pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% adjust nobs to feed the lags
nobse = nobs - nlag;

ylag = glags(y, nlag);
% reslag = glags(resid2,nlag);
% vp=y(nlag+1:end,1:neqs)-resid2(nlag+1:end,1:neqs);
if (nx > 0)
    xlag = ltflag(x, nlag);
else
    xlag = [];
end


vgam = str.vgam; % the parameters for the constant are at the end of vgam
resid = resid2;
sigmar = sigmar2;
phis = str.phis;
thetast = str.thetast;
reslag = glags(resid, nlag);
vp = y(nlag+1:end, 1:neqs) - resid(nlag+1:end, 1:neqs);

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


% form design matrix for the regression.
zr = [];
for i = 1:length(vgam)
    if isnan(vgam(i))
        zr = [zr, xmat(:, i)];
    end
end


[R, p] = chol(sigmar);
L = R';
if (p > 0)
    error('covariance matrix of residuals2 singular in conmedfjac')
end
thzrt = [];
[nzr, mzr] = size(zr);
Phi0 = phis(:, :, 1);
for i = 1:nobse
    thzrt = [thzrt; Phi0 \ zr((i - 1)*neqs+1:i*neqs, :)];
end
Th = [];
for i = 2:nlag + 1
    Th = [Th, thetast(:, :, i)];
end

V = zeros(nobs*neqs, mzr); % partial derivatives
for i = nlag + 1:nobs
    Vp = [];
    for j = i - 1:-1:i - nlag
        Vp = [Vp; V((j - 1)*neqs+1:j*neqs, :)];
    end
    V((i - 1)*neqs+1:i*neqs, :) = -Th * Vp - thzrt((i - nlag - 1)*neqs+1:(i - nlag)*neqs, :);
end
for i = nlag + 1:nobs
    V((i - 1)*neqs+1:i*neqs, :) = L \ V((i - 1)*neqs+1:i*neqs, :);
end

V = -V(nlag*neqs+1:end, :); % negative derivatives for the regression
lres2 = L \ resid(nlag+1:end, :)';
res2t = vec(lres2); % standardized residuals

fjac = -V;
ff = res2t;
g = fjac' * ff;
