
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Multivariate Model for 8 U.S. series, 1980-1, 2012-8
%
% Model is y_t = X_t*alpha + mu_t + epsilon_t, where mu_t follows the
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

titles = ['VU', 'PVU', 'VN', 'PVN', 'PERMITS', 'STARTS', ...
    'mesesvn', 'r', '', 'VUN', 'PVUN', 'VNN', 'PVNN', ...
    'PERMITN', 'STARTN', 'MN', 'AFFORN'];
data = load(fullfile('data', 'VIVIUSA.dat'));
data = data(:, 1:9);
%eliminate NaNs
data(any(isnan(data)'), :) = [];

data1 = load(fullfile('data', 'PVUsa.dat'));

%replace PVU series with seasonally adjusted series
data(:, 3) = data1;

yy = log(data(:, 2:end));
titles = titles(1, 3:end);

%make two years of missing data. Obs. number 349 = Jan. 2009.
yy(349:349+23, 1) = NaN(24, 1);
% yy(2:2+23,4)=NaN(24,1);

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

%parameter optimization (commented)
smname = 'viviusafun';
% %Levenberg-Marquardt
% info.f=smname; info.tr=1;
% info.tolf=1e-4; info.tolx=sqrt(info.tolf); info.maxit=300; info.nu0=.01; info.jac=0;
% info.prt=2;
% tic
% [x,J,ff,g,iter,conf] = marqdt(info,xv,y,pfix,pvar,xf,stordt,conc,n,r);
% toc

%Levenberg-Marquardt results
x = [0.1436; 1.4686; 0.0901; 1.8531; 1.7653; -1.4023; -0.3352; 0.6239; 0.2444; ...
    0.6179; 0.2772; -0.0004; 0.0849; 0.7147; 0.1527; 0.0705; 1.0326; 0.5658; ...
    0.7013; 1.0620; 0.9852; -0.0034; 0.4399];
xx = x0;
xx(pvar) = x; %estimated parameters

models = modstr_viviusa(y, xx, pfix, pvar, xf, stordt, conc, n, r);
Z = models.Z;
T = models.T;
G = models.G;
H = models.H;
W = models.W;
X = models.X;
%
% get residuals and estimator
%
[F, e, g, M, Pevf] = viviusafun(x, y, pfix, pvar, xf, stordt, conc, n, r);
% compute residual diagnostics
Ss = e' * e;
Ff = F' * F;
ne = length(e); %residual sum of squares
% disp('concentrated parameter:')
nr = 0;
nbeta = 0;
conp = Ss / (ne - nr); %estimated sigma square
% disp('square root of concentrated parameter:')
sconp = sqrt(conp);
% compute prediction error variance (finite sample)
Pevf = Pevf * conp;
% disp('standard error (finite sample)')
SPevf = sqrt(diag(Pevf));
%Pevf (innovations variance) is concentrated out
G = (G * sconp) / SPevf(1);
H = (H * sconp) / SPevf(1);
