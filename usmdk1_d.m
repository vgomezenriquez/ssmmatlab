%Example of estimation of a univariate structural model
%Series is car drivers killed or seriously injured in Great Britain from
%January 1969 to December 1984 (Durbin and Koopman, 2012).
%Two explanatory variables are included in the model, the price of oil and
%the number of kilometers driven.
%

clear
data = load(fullfile('data', 'Seatbelt.dat'));
x = [];
y = data(:, 1);
y1 = data(:, 4); %number of kilometers driven
y2 = data(:, 5); %price of oil
yor = y;
tname = 'Seatbelt';
fname = fullfile('results', 'Seatbelt.txt');
lam = 1; %do not take logs
Y = [y1, y2]; %matrix for regression variables
nreg = 2;
npr = 0; %number of forecasts

lag = 36;
cw = 1.96;
freq = 12;
for dr = 0:1
    for ds = 0:1
        c0 = sacspacdif(y, tname, dr, ds, freq, lag, cw);
        pause
    end
end
close all


%define univariate structural model: trend, trigonometric
%seasonality, and irregular component
comp.level = [1, 0.1, NaN];
comp.seas = [2, .1, NaN];
comp.irreg = [1, .1, NaN];
freq = 12;
comp.freq = freq;
bg_year = 1969;
bg_per = 1;
datei = cal(bg_year, bg_per, freq);
comp.datei = datei;

%copy npr in mpr and make npr zero for estimation
if npr > 0
    mpr = npr;
    npr = 0;
else
    mpr = 0;
end

%create structure and put model into state space form
[str, ferror] = suusm(comp, y, Y, npr);
if ferror > 0
    return
end

%estimate model
%
[result, str] = usmestim(y, str);

disp(' ');
disp('******************** Results from estimation ********************');
disp(' ');
result.tv = result.tv';
mprintr(result)
disp(' ')
disp('(Parmaters are standard deviation ratios)')
fprintf(1, '%s %9.4f\n', 'Concentrated parameter:', sqrt(result.sigma2c));
disp('(concentrated parameter is a standard deviation)')
disp(' ');
disp('***** Estimated variances of disturbances *****');
disp(' ');
fprintf(1, '%s %9.4g\n', 'Level    :', result.xvf(1)^2*result.sigma2c);
fprintf(1, '%s %9.4g\n', 'Seasonal :', result.xvf(2)^2*result.sigma2c);
fprintf(1, '%s %9.4g\n', 'Irregular:', result.sigma2c);
disp('press any key to continue')
pause

%estimated and fixed parameters
xvf = result.xvf;
xf = result.xf;
%t-values of varma estimated parameters are in result.tv
%t-values of estimated regression parameters are in result.tvr
%Note that the standard errors are divided by the concentrated parameter
%(sqrt(result.sigma2c))

%create estimated model
[X, Z, G, W, T, H, ins, ii, ferror] = pr2usm(xvf, xf, str);

disp(' ');
disp('***** Estimated Model  *****');
disp(' ');
disp('Matrix T:')
disp(T)
disp('Matrix Z:')
disp(Z)
disp('Matrix G:')
disp(G)
disp('Matrix H:')
disp(H)
disp('More estimation and diagnostic details are in file "Seatbelt.txt"')
disp('in the subdirectory "results"')
disp('press any key to continue')
pause

%compute recursive residuals
[Xt, Pt, g, M, initf, recrs] = scakff(y, X, Z, G, W, T, H, ins, ii);

%residual diagnostics
e = result.e;
F = result.F;
Ss = e' * e;
Ff = F' * F;
ne = length(e); %residual sum of squares
Pevf = result.Pevf; %prediction error variance
% disp('standard error (finite sample)')
SPevf = result.SPevf;
ny = length(y);
pvar = str.pvar;
nr = length(pvar);
X = str.X;
[junk, nbeta] = size(X);
ndrs = ne + nbeta;
lagl = min(36, max([floor(.2*ny), 3 * freq, 10]));
infr = rescomp(e, lagl, nr, Ss, Pevf, SPevf, Ff, ndrs, nbeta);

%plot residual diagnostics
plotres([], [], [], [], [], 1.96, 'residuals', 1, 0, [], 0, [], infr, 1, 1);

%file for output
fid = fopen(fname, 'w');
% fid=1;
%print estimation results
printusmer(fid, datei, tname, yor, y, ny, lam, str, result, nreg, nbeta);

%print residual diagnostics
printres(fid, infr);
%close external file
if fid ~= 1
    fclose(fid);
end
pause
close all


%smoothing
X = str.X;
W = str.W;
npr = mpr;
if ~isempty(X)
    X = X(1:end-npr, :);
end
if ~isempty(W)
    W = W(1:end-npr, :);
end

[Xt, Pt, g, M] = scakfs(y, X, Z, G, W, T, H, ins, ii);
%smoothing can also be done using the following
% [mh,nh]=size(H); C=eye(mh); D=zeros(mh,nh);
% [mb,nb]=size(X); [mw,nw]=size(X); nb=max(nb,nw); U=zeros(mh,nb); mucd=mh;
% [Xt,Pt,g,M]=smoothgen(y,X,Z,G,W,T,H,ins,ii,mucd,U,C,D);

trend = Xt(:, 1) + X * g(end-1:end);
if (mpr > 0)
    %forecast of trend
    trendp = alpr(1, :)' + Xp * g(end-1:end);
else
    pry = [];
    trendp = [];
end
names = char('Original Series', 'Trend');
tsplot([[y; pry'], [trend; trendp]], datei, names);
pause
close all
