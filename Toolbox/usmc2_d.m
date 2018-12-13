%Example of estimation of a univariate structural model
%The time series data analyzed in this example are annual age-adjusted
%melanoma incidences from the Connecticut Tumor Registry (Houghton,
%Flannery, and Viola 1980) for the years 1936 to 1972. The observations
%represent the number of melanoma cases per 100,000 people.
%

clear
% load data.
yy = load(fullfile('data', 'melanoma.dat'));
x = [];
y = vec(yy');
y = y(1:end-3);
yor = y;
tname = 'melanoma';
lam = 1;
if lam == 0
    yl = log(y); %transform series
    y = yl;
end
Y = []; %matrix for regression variables
npr = 10; %number of forecasts

lag = 12;
cw = 1.96;
freq = 1;
ds = 0;
for dr = 0:1
    c0 = sacspacdif(y, tname, dr, ds, freq, lag, cw);
    pause
end
close all

%define univariate structural model: trend, slope, trigonometric
%seasonality, cycle, irregular and autoregressive component
comp.level = [-1, 0, 0];
comp.slope = [-1, 0, 0];
comp.irreg = [1, .1, NaN];
comp.cycle = [1, .1, NaN];
twopi = 2 * pi;
comp.cyclep = [.9, twopi / 10.; NaN, NaN];
comp.cycleb = [twopi / 15., twopi / 5.];
comp.freq = freq;
bg_year = 1936;
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
disp('press any key to continue')
disp(' ')
fprintf(1, '%s %9.4f\n', 'Concentrated parameter:', sqrt(result.sigma2c));
pause

%estimated and fixed parameters
xvf = result.xvf;
xf = result.xf;
%t-values of varma estimated parameters are in result.tv
%t-values of estimated regression parameters are in result. tvr
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
disp('More estimation and diagnostic details are in file "melanoma.txt"')
disp('in the subdirectory "results"')
disp('press any key to continue')
pause

%compute recursive residuals
[Xt, Pt, g, M, initf, recrs] = scakff(y, X, Z, G, W, T, H, ins, ii);

%residual diagnostics
e = result.e;
F = result.F;
Ss = result.Ss;
Ff = result.Ff;
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
fname = fullfile('results', 'melanoma.txt');
fid = fopen(fname, 'w');
% fid=1;
%print estimation results
nreg = 0;
printusmer(fid, datei, tname, yor, y, ny, lam, str, result, nreg, nbeta);

%print residual diagnostics
printres(fid, infr);
%close external file
if fid ~= 1
    fclose(fid);
end

%compute forecasts
if mpr > 0
    %hb, Mb, A and P are in structure result. Here, hb is the vector of
    %regression estimates and Mb is the matrix of mse. A is the
    %estimated state vector, x_{t|t-1}, obtained with the Kalman filter at the
    %end of the sample and P is the matrix of standard errors.
    hb = result.h;
    Mb = result.M;
    A = result.A;
    P = result.P;
    
    npr = mpr;
    [str, ferror] = suusm(comp, y, Y, npr);
    Xp = str.X;
    Wp = str.W;
    if ~isempty(Xp)
        Xp = Xp(end-npr+1:end, :);
    end
    if ~isempty(Wp)
        Wp = Wp(end-npr+1:end, :);
    end
    cw = 1.96;
    m = 1; %number of series
    [pry, mypr, alpr, malpr] = ssmpred(npr, m, A, P, Xp, Z, G, Wp, T, H, hb, Mb);
    spry = zeros(m, npr);
    sconp = sqrt(result.sigma2c);
    for i = 1:npr
        spry(:, i) = sqrt(diag(mypr(:, :, i))) * sconp;
    end
    %obtain forecasts in the original scale using the log-normal
    %distribution
    opry = pry;
    ospry = spry;
    if lam == 0
        for i = 1:npr
            opry(i) = exp(pry(i)+(spry(i)^2)/double(2.));
            ospry(i) = exp(double(2.)*pry(i)+spry(i)^2) * (exp(spry(i)^2) - double(1.));
        end
    end
    
    %plot forecasts
    out.pry = pry;
    out.spry = spry;
    out.opry = opry;
    out.ospry = ospry;
    out.y = y;
    out.yor = yor;
    out.ny = ny;
    out.npr = npr;
    out.cw = cw;
    out.tname = tname;
    out.lam = lam;
    out.s = freq;
    pfctsusm(out);
end

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
%vector g contains the estimates of the vector (delta',beta')'. Thus, the
%vector of regression estimates, hat(beta), is at the end of g.
%smoothing can also be done using the following
% [mh,nh]=size(H); C=eye(mh); D=zeros(mh,nh);
% [mb,nb]=size(X); [mw,nw]=size(X); nb=max(nb,nw); U=zeros(mh,nb); mucd=mh;
% [Xt,Pt,g,M]=smoothgen(y,X,Z,G,W,T,H,ins,ii,mucd,U,C,D);

cycle = Xt(:, 1);
trend = X * g;
%forecast of cycle and trend
cyclep = alpr(1, :)';
trendp = Xp * g;
tsplot([trend; trendp], datei, 'Trend with forecasts');
pause
tsplot([cycle; cyclep], datei, 'Cycle with forecasts');
pause
close all
