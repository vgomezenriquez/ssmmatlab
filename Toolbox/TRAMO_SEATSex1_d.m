%
%Example of a series treated with TRAMO/SEATS. Series is US IPI from 1960-I
%to  2011-III, quarterly series.
%

clear

% Model identified by TRAMO is (1,1,0)(0,1,1)_4 without mean
% outliers detected by TRAMO;    1960-I, 2011-III
%  61 LS    ( 1 1975)
% The same model and outlier are identified by SSMMATLAB. The paramater
% estimates also coincide.
%
sname = 'USIPI';
out = arimaestos(sname);

%
%SEATS residuals
%
[sres] = SEATSres(out);
%
%ols recursive residuals
%
[olsrress, olsrres] = OLSrres(out);
%
%OLS residuals
%
[olsres] = OLSres(out);
%
%smoothed residuals
%
[smtres] = SMTres(out);

s = out.freq;
if out.model.d > 0
    d = out.model.d;
else
    d = 0;
end
if out.model.ds > 0
    ds = out.model.ds;
else
    ds = 0;
end

%
%plot smoothed and ols residuals.There are n smt residuals and n-ndelta ols
%residuals, where n is the length of the series y and ndelta is the degree
%of differencing.
%
[ns, ms] = size(smtres);
ndelta = d + ds * s;
vnames = char('olsres', 'smtres');
bg_year = out.nziyip(2);
bg_per = out.nziyip(3);
res = cal(bg_year, bg_per, s, ndelta+1);
ndatei = cal(res.year, res.period, s);
figure
tsplot([olsres, smtres(ndelta+1:ns)], ndatei, vnames);
pause
%
%plot standardized ols recursive residuals.There are n-ndelta-nbeta of
%those residuals, where n is the length of the series, ndelta is the degree
%of differencing and nbeta is the number of regressors. These residuals are
%uncorrelated and preserve the temporal dimension.
%
nbeta = length(out.model.hb);
ndelta = d + ds * s;
vnames = 'standardized ols recursive residuals';
res = cal(bg_year, bg_per, s, ndelta+nbeta+1);
ndatei = cal(res.year, res.period, s);
figure
tsplot(olsrress, ndatei, vnames);
pause


%
%Signal extraction using the canonical decomposition
%
Ycomp{1} = 'trend'; %assign the outlier to the trend
outa = arimasigex(out, Ycomp);

%plot spectra
if out.gft == 1
    plotspcd(outa)
end

thrc = outa.compcd.ptnum; % trend-cycle numerator
phir = outa.compcd.ptden; % trend-cycle denominator
%   outa.compcd.ptnur;    % number of nonstationary roots in phir
sigma2r = outa.compcd.ptvar; % variance of the trend-cycle innovations (*)
thsc = outa.compcd.stnum; % seasonal numerator
phis = outa.compcd.stden; % seasonal denominator
%   outa.compcd.stnur; % number of nonstationary roots in phis
sigma2s = outa.compcd.stvar; % variance of the seasonal innovations (*)
thtc = outa.compcd.rt; % transitory component (MA term)
sigma2t = outa.compcd.rtvar; % variance of the transitory component innovations (*)
sigma2i = outa.compcd.itvar; % variance of the irregular component (*)
phitc = outa.compcd.phi; % stationary AR trend polynomial
%(*) in units of the series model innovations

disp('trend-cycle numerator:')
disp(thrc)
disp('trend-cycle denominator:')
disp(phir)
disp('variance of the trend-cycle innovations (*)')
disp(sigma2r)
pause
disp('seasonal numerator:')
disp(thsc)
disp('seasonal denominator:')
disp(phis)
disp('variance of the seasonal innovations (*)')
disp(sigma2s)
pause
disp('stationary AR trend polynomial')
disp(phitc)
pause
disp('variance of the irregular component (*)')
disp(sigma2i)
if ~isempty(thtc)
    disp('transitory numerator:')
    disp(thtc)
    phit = 1.;
    disp('transitory denominator:')
    disp(phit)
    disp('variance of the transitory innovations (*)')
    disp(sigma2t)
end
disp('(*) in units of var(A)')


% Models for the components obtained by SEATS (Canonical decomposition),
% assuming that the model has been estimated by TRAMO
% Model estimated by TRAMO
% phi=[-.49555 1]; Th=[-.86020 1]; stand. dev. of resid.: 0.1661D-01
%
%                   MODELS FOR THE COMPONENTS
%
%  TREND-CYCLE NUMERATOR
%    1.00000000000000      -0.120915955672721      -0.968885573465440
%   0.152030382207280
%  TREND-CYCLE DENOMINATOR
%      1.0000    -2.4956     1.9911    -0.4956
%  INNOV. VAR. (*)     0.31544
%
%  SEAS. NUMERATOR
%    1.00000000000000        1.46717885612683        1.02591246661574
%   1.787979416833018E-002
%  SEAS. DENOMINATOR
%      1.0000     1.0000     1.0000     1.0000
%  INNOV. VAR. (*)     0.00187
%
%  IRREGULAR
%  VAR.     0.09674
%
% (*)   IN UNITS OF VAR(A)
% STANDARD DEVI. OF RESID=  0.1661D-01
%
