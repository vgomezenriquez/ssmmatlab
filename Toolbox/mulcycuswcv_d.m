%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Multivariate Model for 11 U.S. series, 1953-4, 2007-9
% Without cumulator Variables. Flow and average variables are treated as
% missing variables.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;

% [data,titles]=xlsread('data\FINAL_CKZ_JAE2.xls');

data = load(fullfile('data', 'FINAL_CKZ_JAE2.dat'));
titlesm = {'IPI', 'Unrate', 'Man', 'Rsal', 'Rfsal'};
titlesq = {'GDP', 'Prod', 'Hours', 'Cons', 'Inv', 'Infl'};

yy = data(:, 2:end);
[m, n] = size(yy);
ym = yy(:, 1:5); %monthly data
yq = yy(:, 6:n); %quarterly data
l = n - 5; %number of quarterly series

r = 1; %number of slope factors

yqf = yq;

%data for the model
y = [ym, yqf];

freqm = 12; % monthly frequency
freqq = 4; % quarterly frequency


%calendar for the whole data
bg_year = 1947;
bg_per = 1;
datei = cal(bg_year, bg_per, freqm);


%initial date for the data
idate = ical(1953, 4, datei);
%calendar for the data in the exercise
idatei = cal(1953, 4, freqm);
%final date for the data
fdate = ical(2007, 9, datei);
y = y(idate:fdate, :);

% lag=24; cw=1.96;
% for i=1:1
%  tname=titlesm(i);
%  for dr=0:1
%   for ds=0:1
%    c0=sacspacdif(ym(:,i),tname,dr,ds,freqm,lag,cw);
%    pause
%   end
%  end
% end
% close all


%initial parameter values
npr = n + r;
np1 = n + 1;
tn = 2 * n;
tnpr = tn + r;
%first parameter in step is concentrated out
fb = ones(n-1, 1) * .1; %parameters in K matrix
stz = ones(n, 1) * .1; %parameters in D_zeta matrix
step = ones(n-1, 1) * .1; %parameters in D_epsilon matrix
ste = .1; %parameters in D_eta matrix
x0 = [fb', stz', step', ste'];
nx = length(x0);
pfix = []; %fixed parameters
pvar = 1:nx; %free parameters
xv = x0(pvar);
xf = x0(pfix);
%indices for standard deviations
stordt = [1:n, np1:tn, tn + 1:tn + r]; %standard deviations in the model
% conc=1;                            %standard deviation concentrated out
%changed to conc=10 on 27-9-2012. Results more similar to CKZ
conc = 10; %standard deviation concentrated out

smname = 'mulcycusfunwcv';
%parameter optimization

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
[x, J, ff, g, iter, conf] = marqdt(info, xv, y, pfix, pvar, xf, stordt, conc, n, r, l, idatei);
toc
xx = x0;
xx(pvar) = x; %estimated parameters


%
% get residuals and estimator
%
[F, e, g, M, Pevf, A, P] = mulcycusfunwcv(x, y, pfix, pvar, xf, stordt, conc, n, r, l, idatei);

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


fname = fullfile('results', 'mulcycuswcv.txt');
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
modelsf = modstr_mulcycuswcv(x, pfix, pvar, xf, stordt, conc, n, r, l);

%smoothing
[Xt, Pt, g, M] = mulcycussmthwcv(x, y, pfix, pvar, xf, stordt, conc, n, r, l, idatei);


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
tp = 1:ny;

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
tp = 1:654;
%plot data against trends
plot(tp, y(:, 1), tp, Xt(:, 1), 'r'), pause
plot(tp, y(:, 2), tp, Xt(:, 2), 'r'), pause
plot(tp, y(:, 3), tp, Xt(:, 3), 'r'), pause
plot(tp, y(:, 4), tp, Xt(:, 4), 'r'), pause
plot(tp, y(:, 5), tp, Xt(:, 5), 'r'), pause
plot(tp, y(:, 6), tp, Xt(:, 6), 'r'), pause
plot(tp, y(:, 7), tp, Xt(:, 7), 'r'), pause
plot(tp, y(:, 8), tp, Xt(:, 8), 'r'), pause
plot(tp, y(:, 9), tp, Xt(:, 9), 'r'), pause
plot(tp, y(:, 10), tp, Xt(:, 10), 'r'), pause
plot(tp, y(:, 11), tp, Xt(:, 11), 'r'), pause
close all


%save variables
save(fullfile('results', 'mulcycuswcv.mat'))
