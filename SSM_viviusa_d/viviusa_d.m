%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Multivariate Model for 8 U.S. series, 1980-1, 2012-8
%
% The model is y_t = X_t*alpha + mu_t + epsilon_t, where mu_t follows the
% model
%
%   mu_{t+1}   = mu_{t}   + K*beta_{t} + eta_{t}
%   beta_{t+1} = beta_{t} +              zeta_{t},
%
%   K = [1  ]      or      K = [ 1   0  ]
%       [b_2]                  [b_2  1  ]
%       [b_3]                  [b_3  c_3]
%        ...                      ...
%       [b_8]                  [b_8  c_8]
%
% Thus, beta_{t} has dimension one or two. Matrix X for the regression
% variables X_{t} is defined in file modstr_viviusa.m, together with the
% other model matrices.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;

% [data,titles]=xlsread('data\VIVIUSA.xls');
% save('data\VIVIUSA.dat','data','-ascii');
titles = ['VU', 'PVU', 'VN', 'PVN', 'PERMITS', 'STARTS', ...
    'mesesvn', 'r', '', 'VUN', 'PVUN', 'VNN', 'PVNN', ...
    'PERMITN', 'STARTN', 'MN', 'AFFORN'];
save(fullfile('data', 'VIVIUSA.txt'), 'titles', '-ascii');
data = load(fullfile('data', 'VIVIUSA.dat'));
data = data(:, 1:9);
%eliminate NaNs
data(any(isnan(data)'), :) = [];

% [data1,titles1]=xlsread('data\PVUsa.xls');
% save('data\PVUsa.dat','data1','-ascii');
data1 = load(fullfile('data', 'PVUsa.dat'));

%replace PVU series with seasonally adjusted series
data(:, 3) = data1;

yy = log(data(:, 2:end));
titles = titles(1, 3:end);

%make two years of missing data. Obs. number 349 = Jan. 2009.
yy(349:349+23, 1) = NaN(24, 1);


freqm = 12; % monthly frequency


%calendar for the whole data
bg_year = 1980;
bg_per = 1;
datei = cal(bg_year, bg_per, freqm);


%initial date for the data
idate = ical(1980, 1, datei);
%final date for the data
fdate = ical(2012, 9, datei);
y = yy(idate:fdate, :);
[m, n] = size(y);


% analyze individual series
% lag=24; cw=1.96;
% for i=1:n
%  tname=titles(i);
%  for dr=0:1
%   for ds=0:1
%    c0=sacspacdif(y(:,i),tname,dr,ds,freqm,lag,cw);
%    pause
%   end
%  end
% end
% close all


r = 1;
idatei = datei;
%initial parameter values
np1 = n + 1;
tn = 2 * n;
tnpr = tn + r;
%number of parameters in K matrix
nk = 0;
for i = 1:r
    nk = nk + n - i;
end
fb = ones(nk, 1) * .1; %parameters in K matrix
%one standard deviation is concentrated out
stz = ones(n-1, 1) * .1; %parameters in D_zeta matrix
step = ones(n, 1) * .1; %parameters in D_epsilon matrix
ste = ones(r, 1) * .1; %parameters in D_eta matrix
x0 = [fb', stz', step', ste'];
nx = length(x0);
pfix = []; %fixed parameters
pvar = 1:nx; %free parameters
xv = x0(pvar);
xf = x0(pfix);
%indices for standard deviations
stordt = [1:n, np1:tn, tn + 1:tnpr]; %standard deviations in the model
conc = 7; %standard deviation concentrated out


%parameter optimization
smname = 'viviusafun';

%Levenberg-Marquardt
info.f = smname;
info.tr = 1;
info.tolf = 1e-4;
info.tolx = sqrt(info.tolf);
info.maxit = 300;
info.nu0 = .01;
info.jac = 0;
info.prt = 2;
tic
[x, J, ff, g, iter, conf] = marqdt(info, xv, y, pfix, pvar, xf, stordt, conc, n, r);
toc
xx = x0;
xx(pvar) = x; %estimated parameters


%
% get residuals and estimator
%
[F, e, g, M, Pevf, A, P] = viviusafun(x, y, pfix, pvar, xf, stordt, conc, n, r);

% compute residual diagnostics
Ss = e' * e;
Ff = F' * F;
ne = length(e); %residual sum of squares
disp('concentrated parameter:')
nr = 0;
nbeta = 0;
conp = Ss / (ne - nr); %estimated sigma square
disp('square root of concentrated parameter:')
sconp = sqrt(conp);
% compute prediction error variance (finite sample)
disp('prediction error variance (finite sample)')
Pevf = Pevf * conp;
disp('standard error (finite sample)')
SPevf = sqrt(diag(Pevf));
%
% lagl=min(36,max([floor(.2*ny) 3*s 10]));
lagl = 36;
ndrs = ne + nbeta;
infr = rescomp(e, lagl, nr, Ss, conp, sconp, Ff, ndrs, nbeta);


% %standard errors via second derivatives of log likelihood
% xt=x';
% H=fdhess('logF',xt,smname,y,pfix,pvar,xf,stordt,conc,n,r,l,idatei);
% SS=inv(H/2)/(ne-nr); se=sqrt(abs(diag(SS)))';
% %t-values
% disp('t-values:')
% tp=zeros(size(xx));
% tt(pfix)=NaN;
% tt(pvar)=x./se

%standard errors and t-values
V = -J; % negative derivatives for the regression
beta = x';
res2t = ff + V * beta; % regression
[beta3, tv] = btval([], [res2t, V]);
tt(pfix) = NaN;
tt(pvar) = tv


fname = fullfile('results', 'viviusa.txt');
%file for output
fid = fopen(fname, 'w');

%print residual diagnostics
printres(fid, infr);
%end of residual diagnostics

if fid ~= 1
    fclose(fid);
end


%plot residual diagnostics
plotres([], [], [], [], [], 1.96, 'residuals', 1, 0, [], 0, [], infr, 1, 1);
%close figures
close all

%estimated model
modelsf = modstr_viviusa(y, x, pfix, pvar, xf, stordt, conc, n, r);

%smoothing
[Xt, Pt, g, M] = viviusasmth(x, y, pfix, pvar, xf, stordt, conc, n, r);


%Mse for the trends
[ny, my] = size(y);
[nx, mx] = size(Xt);
strend = zeros(size(y));
for i = 1:ny
    ia = (i - 1) * mx + 1:i * mx;
    P = Pt(ia, 1:mx);
    for j = 1:my
        strend(i, j) = sconp * sqrt(P(j, j));
    end
end
% confidence interval width: cw=1.69
cw = 1.69;

bg_yr1 = 1953;
bg_per1 = 4;
datei1 = cal(bg_yr1, bg_per1, freqm);

%all trends with 95% intervals and dates
for i = 1:my
    figure
    trend = Xt(:, i);
    trendpcw = trend + cw * strend(:, i);
    trendmcw = trend - cw * strend(:, i);
    vnames = char(['trend', num2str(i), '+.95mse'], ['series', num2str(i)], ...
        ['trend', num2str(i)], ['trend', num2str(i), '-.95mse']);
    tsplot([trendpcw, y(:, i), trend, trendmcw], datei1, vnames), pause
end
close all


%data length
tp = 1:m;
%plot data against trends
for i = 1:n
    plot(tp, y(:, i), tp, Xt(:, i), 'r'), pause
end
close all

%Forecasting

%model matrices
X = modelsf.X;
Z = modelsf.Z;
G = modelsf.G;
W = modelsf.W;
T = modelsf.T;
H = modelsf.H;
ins = modelsf.ins;
ii = modelsf.i;


%extend series with one year NaNs for forecasting
npr = 16;
yf = NaN(npr, my);
[ny, p] = size(y);
y = [y; yf];

% %extend regression variables for variable VU according to the results
% %given by TRAMO. The interventions are: TC 367, LS 360, AO 37.
% xnp=npr*n;
% Xx1=zeros(xnp,3); fac=X((ny-1)*p+1,1);
% for j=1:npr
%  ip=(j-1)*p+1;
% %intervention TC 367
%   Xx1(ip,1)=fac;
%   fac=min(1.d-10,fac*fac);
% %intervention LS 360
%   Xx1(ip,2)=1.;
% end
Xx1 = [];

%extend regression variables for variable PVN according to the results
%given by TRAMO. The interventions are: AO 370, LS 89, LS 344, LS 127, LS
%289.
xnp = npr * n;
Xx2 = zeros(xnp, 5);
for j = 1:npr
    ip = (j - 1) * p + 4;
    %interventions: LS 89, LS 344, LS 127, LS 289.
    Xx2(ip, 2:5) = ones(1, 4);
end

Xx = [Xx1, Xx2];
[junk, nxpr] = size(Xx);
X = [X; Xx];


fdate = ical(2014, 1, datei);


% %square root filter for smoothing.
% [Xsct,Psct,gsc,Msc]=scakfssqrt(y,X,Z,G,W,T,H,ins,ii);

%forecasts
[Xsct, Psct, gsc, Msc] = smoothgen(y, X, Z, G, W, T, H, ins, ii, n, [], Z, G);

%save forecasts in a spreadsheet
% xlswrite(fullfile('results','viviusaf.xls'),Xsct);

%calendar for the plots
bg_yr1 = 1980;
bg_per1 = 1;
datei1 = cal(bg_yr1, bg_per1, freqm);
%date format
ydigit = 'yyyy';


% %all trends plus forecasts, whole sample.
% for i=1:my
% figure
% trend=Xsct(:,i); ser=y(:,i);
% vnames=vertcat(titles(i),strcat(titles(i), ' + 12 Forecasts'));
% hold on
% tsplot([ser trend],datei1,vnames,ydigit)
% hold on
% ylim = get(gca,'YLim');
% plot(ones(2,1)*datenum(2012,8,1),ylim,':k')
% hold off
% pause
% end
% close all

%all trends plus forecasts, origin of forecast 2012:8
iobs = ical(2011, 6, datei);
fobs = ical(2014, 1, datei);
bg_yr1 = 2011;
bg_per1 = 6;
datei2 = cal(bg_yr1, bg_per1, freqm);
for i = 1:my
    figure
    trend = Xsct(iobs:fobs, i);
    ser = y(iobs:fobs, i);
    %uncomment next line if there are regression variables
    trend = trend + X((iobs - 1)*p+i:p:(fobs - 1)*p+i, :) * gsc(end-nxpr+1:end);
    vnames = char(titles(i), strcat(titles(i), ' + 16 Forecasts'));
    hold on
    tsplot([ser, trend], datei2, vnames, ydigit)
    hold on
    ylim = get(gca, 'YLim');
    plot(ones(2, 1)*datenum(2012, 9, 1), ylim, ':k')
    hold off
    pause
end
close all
