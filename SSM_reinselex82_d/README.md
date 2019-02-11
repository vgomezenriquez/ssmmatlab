[<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/banner.png" width="888" alt="Visit QuantNet">](http://quantlet.de/)

## [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/qloqo.png" alt="Visit QuantNet">](http://quantlet.de/) **SSM_reinselex82_d** [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/QN2.png" width="60" alt="Visit QuantNet 2.0">](http://quantlet.de/)

```yaml

Name of QuantLet: SSM_reinselex82_d

Published in: Linear Time Series With MATLAB and Octave

Description: 'A VAR and a VARMA model for the data in Exercise 8.2 of Reinsel (1997) are estimated'

Keywords: time-series, VAR model, VARMA model, estimation, stepwise elimination

Author: Víctor Gómez

Submitted: Fri, December 28 2018 by Víctor Gómez

```

### MATLAB Code
```matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example 8.2 of Reinsel (1997), pp. 292-298
%
% Series are: 1) the in-phase current, 2) the out-of-phase current and 3)
% the frequency of the voltage generated, expressed as deviations from
% nominal values. The number of observations is n=100. The data have been
% multiplied by 10.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;

y = load(fullfile('data', 'power-turbo.dat'));

lag = 20;
cw = 1.96;
freq = 1;
ds = 0;
tname = {'In-phase current deviations', 'Out-of-phase current deviations', ...
    'Frequency of voltage deviations'};
for i = 1:3
    for dr = 0:1
        c0 = sacspacdif(y(:, i), tname(i), dr, ds, freq, lag, cw);
        pause
    end
end
close all


%Preliminary VAR analysis
prt = 1;
minlag = 0;
maxlag = 6;
lagsoptlr = varident(y, maxlag, minlag, prt);
disp('Preliminary VAR analysis: lag lenght estimated with LR')
disp(lagsoptlr)
pause

%We estimate a VAR of order 4 as in Reinsel p.295
disp('estimation of a VAR of order 4:')
disp('press any key to continue')
pause
test = 1;
lags = 4;
res = var_est(y, lags, test);

disp(' ');
disp('***** Estimated VAR Model  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'AR';
strt = 1;
mprintar(res.phi(:, :, 2:5), in, tit, strt);
disp(' ')
tit = 'Constant';
mprintar(res.const', in, tit);
disp(' ')
tit = 'Sigma:';
mprintar(res.sigmar, in, tit);

disp(' ');
disp('***** Estimated t-values  *****');
disp(' ');
tit = 'tv-AR';
strt = 1;
mprintar(res.phitv(:, :, 2:5), in, tit, strt);
disp(' ')
tit = 'tv-Constant';
mprintar(res.consttv', in, tit);
disp(' ')
tit = 'Granger caus. prob.:';
mprintar(res.fprob, in, tit);

disp(' ')
disp('press any key to continue')
pause

disp(' ')
disp('estimation of a simplified VAR of order 4 by the HR method ')
disp('after eliminating all (1,3) and (2,3) elements of the matrices')
disp(' Phi_j  as in Reinsel (1997):')
disp('press any key to continue')
pause

seas = 1;
x = [];
hr3 = 1;
finv2 = 1;
[strv, ferror] = estvarmaxpqrPQR(y, x, seas, [4, 0, 0], [0, 0, 0], hr3, finv2);

for i = 2:5
    strv.phi(1, 3, i) = 0;
    strv.nparm = strv.nparm - 1;
    strv.phi(2, 3, i) = 0;
    strv.nparm = strv.nparm - 1;
end
strv = mhanris(y, x, seas, strv, hr3, finv2);

disp(' ');
disp('***** Estimated simplified VAR(4) Model  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'AR';
strt = 1;
mprintar(strv.phis(:, :, 2:5), in, tit, strt);
disp(' ')
tit = 'Constant';
mprintar(strv.mu', in, tit);
disp(' ')
tit = 'tv-AR';
strt = 1;
mprintar(strv.phitv(:, :, 2:5), in, tit, strt);
disp(' ')
tit = 'tv-Constant';
mprintar(strv.mutv', in, tit);
disp(' ')
tit = 'Sigma';
mprintar(strv.sigmar2, in, tit);
disp('press any key to continue')
pause

disp(' ')
disp(' p-value to test exogeneity:')
%statistic to test the constraints
[T, k] = size(y);
m = 4;
qstat = -((T - m) - m * k - 1 - .5) * log(det(res.sigmar)/det(strv.sigmar2));
df = 8;
pval = 1 - gammp(df*.5, qstat*.5)
disp('press any key to continue')
pause

recrs = strv.resid2;
lag = 12;
ic = 1;
disp(' ')
disp('******** VAR Residuals:     ********');
str = mautcov(recrs, lag, ic);
disp('Correlation matrix at lag 0:')
disp(str.r0)
disp('press any key to continue')
pause


disp(' ')
disp('identify a VARMA(p,q) model for the series')
disp(' ')
disp('press any key to continue')
pause
%identify a VARMA(p,q) model for the series
maxlag = 6;
minlag = 0;
prt = 0;
x = [];
seas = 1;
[lagsopt, ferror] = lratiopqr(y, x, seas, maxlag, minlag, prt);
disp(' ')
disp('Estimated orders in VARMAX(p,q,r):  ')
disp(lagsopt)
disp('press any key to continue')
pause


disp(' ')
disp('estimate a constrained VARMA(4,1) using the HR method')
disp('press any key to continue')
pause

%estimate a simplified VARMA(4,1) by the Hannan-Rissanen method
clear strv
seas = 1;
x = [];
hr3 = 0;
finv2 = 1;
[strv, ferror] = estvarmaxpqrPQR(y, x, seas, [4, 1, 0], [0, 0, 0], hr3, finv2);

for i = 2:5
    strv.phi(1, 3, i) = 0;
    strv.nparm = strv.nparm - 1;
    strv.phi(2, 3, i) = 0;
    strv.nparm = strv.nparm - 1;
end
strv = mhanris(y, x, seas, strv, hr3, finv2);


disp(' ');
disp('***** Estimated constrained VARMA(4,1) Model  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'phi';
strt = 1;
mprintar(strv.phis3(:, :, 2:5), in, tit, strt);
disp(' ')
tit = 'th';
strt = 1;
mprintar(strv.thetas3(:, :, 2), in, tit, strt);
disp(' ')
tit = 'Constant';
mprintar(strv.mus3', in, tit);
disp(' ')
tit = 'tv-phi';
strt = 1;
mprintar(strv.phitv3(:, :, 2:5), in, tit, strt);
disp(' ')
tit = 'tv-th';
strt = 1;
mprintar(strv.thetatv3(:, :, 2), in, tit, strt);
disp(' ')
tit = 'tv-Constant';
mprintar(strv.mutv3', in, tit);
disp(' ')
tit = 'Sigma';
mprintar(strv.sigmar3, in, tit);
disp('press any key to continue')
pause

disp(' ')
disp('impose additional constraints in the  VARMA(4,1) model ')
disp('as in Reinsel (1997)  and estimate again using the HR method')
disp('press any key to continue')
pause

%impose additional constraints as in Reinsel (1997)
strv.phi(2, 1, 2) = 0;
strv.phi(3, 2, 2) = 0;
strv.phi(2, 1, 3) = 0;
strv.phi(3, 2, 3) = 0;
strv.phi(3, 3, 3) = 0;
strv.phi(3, 3, 4) = 0;
strv.phi(3, 2, 5) = 0;
strv.theta(3, :, 2) = zeros(1, 3);
strv.theta(2, 1, 2) = 0;
strv.theta(1, 2, 2) = 0;
strv.theta(1:2, 3, 2) = zeros(2, 1);
strv.nparm = strv.nparm - 14;
strv = mhanris(y, x, seas, strv, hr3, finv2);
disp(' ');
disp('***** Estimated constrained VARMA(4,1) Model  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'phi';
strt = 1;
mprintar(strv.phis3(:, :, 2:5), in, tit, strt);
disp(' ')
tit = 'th';
strt = 1;
mprintar(strv.thetas3(:, :, 2), in, tit, strt);
disp(' ')
tit = 'Constant';
mprintar(strv.mus3', in, tit);
disp(' ')
tit = 'tv-phi';
strt = 1;
mprintar(strv.phitv3(:, :, 2:5), in, tit, strt);
disp(' ')
tit = 'tv-th';
strt = 1;
mprintar(strv.thetatv3(:, :, 2), in, tit, strt);
disp(' ')
tit = 'tv-Constant';
mprintar(strv.mutv3', in, tit);
disp(' ')
tit = 'Sigma';
mprintar(strv.sigmar3, in, tit);
disp('press any key to continue')
pause


%estimate using exact ML
%setup model
Phi = eye(3);
Th = eye(3);
phi = strv.phis3;
th = strv.thetas3(:, :, 1:2);
Sigma = strv.sigmar3;
% phi=strc.phiscon; th=strc.thetascon(:,:,1:2); Sigma=strc.sigmarcon;


%create regression variable for the mean
Y = eye(3);
freq = 1;
[str, ferror] = suvarmapqPQ(phi, th, Phi, Th, Sigma, freq);
str.phin = strv.phi;
str.thn = strv.theta(:, :, 1:2);
[str, ferror] = fixvarmapqPQ(str);


disp(' ')
disp('estimation using the exact method')
disp('press any key to continue')
pause

%estimate model using the exact method
result = varmapqPQestim(y, str, Y);
disp(' ');
disp('******************** Results from estimation ********************');
disp(' ');
mprintr(result)
disp('press any key to continue')
pause


%estimated and fixed parameters
xvf = result.xvf;
xf = result.xf;
%t-values of varma estimated parameters are in result.tv
%t-values of estimated regression parameters are in result. tvr

%create estimated model
[phif, thf, Phif, Thf, Lf, ferror] = pr2varmapqPQ(xvf, xf, str);
Sigmar = result.Sigmar;
%t-values
tvf = result.tv;
[phitvf, thtvf, Phitvf, Thtvf, Ltvf, ferror] = pr2varmapqPQ(tvf, xf, str);

disp(' ');
disp('***** Estimated Model using the exact method  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'mu:';
mprintar(result.h', in, tit);

tit = 'phi';
strt = 1;
mprintar(phif(:, :, 2:5), in, tit, strt);

tit = 'th';
strt = 1;
mprintar(thf(:, :, 2), in, tit, strt);
disp(' ');
tit = 'Sigma:';
mprintar(Sigmar, in, tit);
disp(' ')
tit = 'L:';
mprintar(Lf, in, tit);
disp(' ')


disp('t-values: ')
disp(' ');
tit = 'tv-mu:';
mprintar(result.tvr', in, tit);

tit = 'tv-phi';
strt = 1;
mprintar(phitvf(:, :, 2:5), in, tit, strt);

tit = 'tv-th';
strt = 1;
mprintar(thtvf(:, :, 2), in, tit, strt);
tit = 'tv-L:';
mprintar(Ltvf, in, tit);
disp(' ')

disp('press any key to continue')
pause


%compute exact OLS residuals using fast square root filter
Sigmaf = result.Sigmar;
[strf, ferror] = suvarmapqPQ(phif, thf, Phif, Thf, Sigmaf, freq);
s = size(y, 2);
pr = 4;
kro = pr * ones(1, s);
m = 0;
[strr, ferror] = matechelon(kro, s, m);
strr.phis = phif;
strr.thetas = thf;
strr.gammas = [];
strr.thetas(:, :, pr+1) = zeros(s);
strr = armaxe2sse(strr);
strr.sigmar2 = strf.Sigma;

Y = eye(s);
tol = 1.d-10;
maxupdt = [];
[e, E, rSigmat] = sqrt_ckms(y, Y, strr, maxupdt, tol);
[ne, me] = size(e);
recr = zeros(ne, me);
nbeta = s;
for ii = 1:ne
    ind = (ii - 1) * nbeta + 1:ii * nbeta; % V=rSigmat(ind,:);
    recr(ii, :) = e(ii, :) - (E(ind, :)' * result.h)';
end

% %compute OLS residuals using square root filter
% %set up regression matrices
% X=eye(s); W=[];
% Sigmax=strr.sigmar2;
% [L,p] = chol(Sigmax,'lower'); Lm1=pinv(L);
% %set up system matrices
% T=strr.Fs; Z=strr.Hs; G=Lm1; H=strr.Ks*Lm1;
% %set up initial conditions
% ndelta=0;                            %number of unit roots
% [ins,i,ferror]=incossm(T,H,ndelta);
% [KKP,Pt,recrs,recr]=scakfffsqrt(y,X,Z,G,W,T,H,ins,i,g);

%plot residuals
plot(recr(:, 1)), legend('recr(:,1)'), pause
plot(recr(:, 2)), legend('recr(:,2)'), pause
plot(recr(:, 3)), legend('recr(:,3)'), pause
close all
%compute autocovariance and autocorrelation matrices of residuals
lag = 8;
ic = 1;
nr = 0;
disp(' ')
disp('******** OLS Residuals:     ********');
stre = mautcov(recr, lag, ic, nr);
disp('Correlation matrix at lag 0:')
disp(stre.r0)
disp('Q statistics:')
disp(stre.qstat)

disp('p-values of Q statistics:')
disp(stre.pval)
[m, n] = size(stre.pval);
t = 1:m;
plot(t, stre.pval, t, 0.05*ones(1, m))
legend('p-values of Q statistics:')
disp('press any key to continue')
pause
close all


disp(' ')
disp('estimation of a VARMAX(4,0,4) by the Hannan-Rissanen method')
disp('press any key to continue')
pause

%estimate a VARMAX(4,0,4) by the Hannan-Rissanen method. Var y_t3 is the
%output and y_1t and y_2t are the inputs. See p. 297 in Reinsel (1997)
x = y(:, 1:2);
yo = y(:, 3);
[strv, ferror] = estvarmaxpqrPQR(yo, x, seas, [4, 0, 4], [0, 0, 0], hr3, finv2);
strv.phi(1, 1, 3) = 0;
strv.phi(1, 1, 4) = 0;
strv.gamma(1, 1, 1) = 0;
strv.gamma(1, 2, 1) = 0;
strv.gamma(1, 2, 2) = 0;
strv.gamma(1, 2, 3) = 0;
strv.gamma(1, 2, 5) = 0;
strv.nparm = strv.nparm - 7;
strv = mhanris(yo, x, seas, strv, hr3, finv2);

disp(' ');
disp('***** Estimated constrained VARMA(4,1) Model  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'phi';
strt = 1;
mprintar(strv.phis3(:, :, 2:5), in, tit, strt);
disp(' ')
tit = 'gamma';
strt = 0;
mprintar(strv.gammas3(:, :, 1:5), in, tit, strt);
disp(' ')
tit = 'Constant';
mprintar(strv.mus3', in, tit);
disp(' ')
tit = 'tv-phi';
strt = 1;
mprintar(strv.phitv3(:, :, 2:5), in, tit, strt);
disp(' ')
tit = 'tv-gamma';
strt = 0;
mprintar(strv.gammatv3(:, :, 1:5), in, tit, strt);
disp(' ')
tit = 'tv-Constant';
mprintar(strv.mutv3', in, tit);
disp(' ')
tit = 'Sigma';
mprintar(strv.sigmar3, in, tit);
disp('press any key to continue')
pause


disp(' ')
disp('estimation using the exact method')
disp('press any key to continue')
pause

Y = 1;
[xvf, strx, ferror] = mexactestim(yo, x, strv, Y);


disp(' ');
disp('***** Estimated Model using the exact method  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'mu:';
mprintar(strx.musexct', in, tit);

tit = 'phi';
strt = 1;
mprintar(strx.phisexct(:, :, 2:5), in, tit, strt);

tit = 'gamma';
strt = 0;
mprintar(strx.gammasexct(:, :, 1:5), in, tit, strt);
disp(' ');
tit = 'Sigma:';
mprintar(strx.sigmarexct, in, tit);
disp(' ')


disp('t-values: ')
disp(' ');
tit = 'tv-mu:';
mprintar(strx.mutvexct', in, tit);

tit = 'tv-phi';
strt = 1;
mprintar(strx.phitvexct(:, :, 2:5), in, tit, strt);

tit = 'tv-gamma';
strt = 0;
mprintar(strx.gammatvexct(:, :, 1:5), in, tit, strt);
disp('press any key to continue')
pause

disp(' ')
disp('check the number of unit roots in the model')
disp(' ')
[D, nr, yd, DA, ferror] = mcrcregr(y);
disp(' ')
disp('number of unit roots according to the crc criterion:')
disp(nr)

```

automatically created on 2019-02-11