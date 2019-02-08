%Example of estimation of an ARMA(p,q)(P,Q)_s model
%Series is airline series from Box and Jenkins (1976)
%

clear
% load data.
y = load(fullfile('data', 'bjsgairl.dat'));
x = [];
yl = log(y); %transform series
Y = []; %matrix for regression variables


lag = 36;
cw = 1.96;
freq = 12;
tname = {'BJ Airline Passengers'};
for dr = 0:1
    for ds = 0:1
        c0 = sacspacdif(yl, tname, dr, ds, freq, lag, cw);
        pause
    end
end
close all

%model identification

%1) select differencing degrees using CRC criterion
maxr = 2;
[nr1, ns1, nr, ns] = crcreg(yl, freq, maxr);

dr = nr;
ds = ns;
yd = diffest(yl, Y, freq, 0, dr, ds, 0, 0); %differenced series

%2) select ARMA degrees using BIC
parm.s = freq;
parm.S = 0;
%intial model
parm.dr = dr;
parm.ds = ds;
parm.dS = 0;
parm.p = 0;
parm.ps = 0;
parm.q = 0;
parm.qs = 0;
parm.qS = 0;
%next two paramaters for Hannan-Rissanen method
ols = 0;
a = 2.;
maxpq = 3; %maximum p,  q degrees
maxPQ = 1; %maximum ps, qs degrees

parm = armaid(yl, parm, ols, a, maxpq, maxPQ);
%model selected
p = parm.p;
ps = parm.ps;
q = parm.q;
qs = parm.qs;
%end of model identification

disp(' ');
disp('******************** Model selected ********************');
disp(' ');
parm
disp('press any key to continue')
pause

% define model. Model is (0,0,1)(0,0,1)_12 for the differenced logged
% series
phi(:, :, 1) = 1;
Phi(:, :, 1) = 1;
th(:, :, 1) = 1;
Th(:, :, 1) = 1;
%no mean in the model
Y = [];
npr = 12; %number of forecasts
%copy npr in mpr and make npr zero for estimation
if npr > 0
    mpr = npr;
    npr = 0;
else
    mpr = 0;
end

%estimate model using HR method
[strv, ferror] = estvarmaxpqrPQR(yd, x, freq, [p, q, 0], [ps, qs, 0], 0, 1, 1);

%setup model
th(:, :, 2) = strv.thetas3(:, :, 2);
Th(:, :, 2) = strv.thetas3(:, :, freq+1);
Sigma = strv.sigmar3;


%create structure and put model into state space form
[str, ferror] = suvarmapqPQ(phi, th, Phi, Th, Sigma, freq);

%estimate model
result = varmapqPQestim(yd, str, Y);

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

disp(' ');
disp('***** Estimated Model  *****');
disp(' ');
clear in
in.fid = 1;
in.cnames = char(' th(1):', ' Th(1):', ' Sigma:');
in.fmt = char('%12.4f');
z = [thf(:, :, 2), Thf(:, :, 2), result.Sigmar];
mprint(z, in);
disp('press any key to continue')
pause


%residual diagnostics
e = result.e; %white noise residuals
ff = result.ff; %vector of nonlinear functions
nbeta = 0; %length of regression vector
ndrs = length(yd);
Ss = e' * e; %residual sum of squares
Ff = ff' * ff;
conp = result.sigma2c;
sconp = sqrt(conp); %residual standar error
lagl = 3 * freq;
infr = rescomp(e, lagl, length(xvf), Ss, conp, sconp, Ff, ndrs, nbeta);

ic = 1;
lag = 26;
disp(' ')
disp('******** Residuals:     ********');
str = mautcov(e, lag, ic);
disp('More estimation and diagnostic details are in file "bjsgairl.txt"')
disp('in the subdirectory "results"')
disp('press any key to continue')
pause


%plot residual diagnostics
plotres([], [], [], [], [], 1.96, 'residuals', 1, 0, [], 0, [], infr, 1, 1);
close all

%print residual diagnostics
%file for output
fname = fullfile('results', 'bjsgairl.txt');
fid = fopen(fname, 'w');
% fid=1;
printres(fid, infr);
%close external file
if fid ~= 1
    fclose(fid);
end

%compute forecasts of the logged differenced series
if mpr > 0
    %hb, Mb, A and P are in structure result. Here, hb is the vector of
    %regression estimates and Mb is the matrix of standard errors. A is the
    %estimated state vector, x_{t|t-1}, obtained with the Kalman filter at the
    %end of the sample and P is the matrix of standard errors.
    hb = result.h;
    Mb = result.H;
    A = result.A;
    P = result.P;
    
    npr = mpr;
    %set up system matrices for the estimated model
    %Note that the residual covariance matrix is divided by the concentrated
    %parameter (result.sigma2c).
    Sigmaf = Lf * Lf';
    [strf, ferror] = suvarmapqPQ(phif, thf, Phif, Thf, Sigmaf, freq);
    Z = strf.Z;
    G = strf.G;
    T = strf.T;
    H = strf.H;
    Xp = Y;
    Wp = [];
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
    
    %plot forecasts
    tname = 'bjsgairl (differenced and in logs)';
    out.pry = pry;
    out.spry = spry;
    out.opry = opry;
    out.ospry = ospry;
    out.y = yd;
    out.yor = yd;
    out.ny = length(yd);
    out.npr = npr;
    out.cw = cw;
    out.tname = tname;
    lam = 1; %lam=0, logs are taken; =1, no logs are taken
    %in this case, since we work with the logged
    %series, lam=1.
    out.lam = lam;
    out.s = freq;
    pfctsusm(out);
    
end

%compute forecasts of the original series
if mpr > 0
    npr = mpr;
    %set up system matrices for the estimated ARIMA model
    %Note that the residual covariance matrix is divided by the concentrated
    %parameter (result.sigma2c).
    Sigmaf = Lf * Lf';
    %Differencing polynomial
    phifo(:, :, 1) = 1.;
    phifo(:, :, 2) = -1.;
    Phifo(:, :, 1) = 1.;
    Phifo(:, :, 2) = -1.;
    %MA polynomial
    thfo = thf;
    Thfo = Thf;
    [strfo, ferror] = suvarmapqPQ(phifo, thfo, Phifo, Thfo, Sigmaf, freq);
    %ARIMA model in state space form
    Z = strfo.Z;
    G = strfo.G;
    T = strfo.T;
    H = strfo.H;
    [ndelta, junk] = size(T);
    X = [];
    W = [];
    %initial conditions for the Kalman filter
    [ins, i, ferror] = incossm(T, H, ndelta);
    chb = 0; %there are no regression effects, so do not compute hb and Mb in
    %scakfle2
    
    %run Kalman filter
    [e, f, hb, Mb, A, P, qyy, R] = scakfle2(yl, X, Z, G, W, T, H, ins, i, chb);
    %hb is the vector of regression estimates and Mb is the matrix of standard
    %errors. A is the estimated state vector, x_{t|t-1}, obtained with the
    %Kalman filter at the end of the sample and P is the matrix of standard
    %errors.
    
    %forecasts
    [pry, mypr, alpr, malpr] = ssmpred(npr, m, A, P, Xp, Z, G, Wp, T, H, hb, Mb);
    spry = zeros(m, npr);
    sconp = sqrt(result.sigma2c);
    for i = 1:npr
        spry(:, i) = sqrt(diag(mypr(:, :, i))) * sconp;
    end
    %obtain forecasts in the original scale using the log-normal
    %distribution
    lam = 0;
    opry = pry;
    ospry = spry;
    if lam == 0
        for i = 1:npr
            opry(i) = exp(pry(i)+(spry(i)^2)/double(2.));
            ospry(i) = exp(double(2.)*pry(i)+spry(i)^2) * (exp(spry(i)^2) - double(1.));
        end
    end
    %plot forecasts
    tname = 'bjsgairl';
    out.pry = pry;
    out.spry = spry;
    out.opry = opry;
    out.ospry = ospry;
    out.y = yl;
    out.yor = y;
    out.ny = length(yl);
    out.npr = npr;
    out.cw = cw;
    out.tname = tname;
    out.lam = lam;
    out.s = freq;
    pfctsusm(out);
end
