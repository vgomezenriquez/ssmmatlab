% Example of an ARIMA(p,d,q) model
% Series is series e from Box and Jenkins (1976)
% WOLFER SUNSPOT NUMBERS: YEARLY, 1770 - 1869
%

clear
y = load(fullfile('data', 'seriee.dat'));
x = [];
yor = y;
tname = 'seriee';
p = 2;
d = 0;
q = 1;
lam = 1;
if lam == 0
    yl = log(y); %transform series
    y = yl;
end
n = size(y, 1);
npr = 10; %number of forecasts
Y = ones(n+npr, 1); %matrix for regression variables


%copy npr in mpr and make npr zero for estimation
if npr > 0
    mpr = npr;
    npr = 0;
    YY = Y;
    Y = Y(1:n, :);
else
    mpr = 0;
end

lag = 20;
cw = 1.96;
freq = 1;
ds = 0;
for dr = 0:d
    c0 = sacspacdif(y, tname, dr, ds, freq, lag, cw);
    pause
end
close all


%
%we specify an ARIMA(p,d,q) model
yd = y;
for i = 1:d
    yd = diferm(yd, 1); %differenced series
end

%estimate model using HR method
[strv, ferror] = estvarmaxpqrPQR(yd, x, freq, [p, q, 0], [0, 0, 0], 0, 1, 1);


% define ARMA model for the differenced series. Model is (p,q)
phi(:, :, 1) = 1;
Phi(:, :, 1) = 1;
th(:, :, 1) = 1;
Th(:, :, 1) = 1;

%setup model for the differenced series
for i = 1:p
    phi(:, :, i+1) = strv.phis3(:, :, i+1);
end
for i = 1:q
    th(:, :, i+1) = strv.thetas3(:, :, i+1);
end
Sigma = strv.sigmar3;


%create structure for the differenced series and put model into state space
%form
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
Sigmar = result.Sigmar;

disp(' ');
disp('***** Estimated Model  *****');
disp(' ');
clear in
in.fid = 1;
if (p > 0)
    in.cnames = ' phi(1):';
    for i = 2:p
        in.cnames = char(in.cnames, [' phi(', num2str(i), '):']);
    end
end
if (q > 0)
    if (p > 0)
        in.cnames = char(in.cnames, ' th(1):');
    else
        in.cnames = ' th(1):';
    end
    for i = 2:q
        in.cnames = char(in.cnames, [' th(', num2str(i), '):']);
    end
end
in.cnames = char(in.cnames, ' Sigma:');
in.fmt = char('%12.4f');
z = [];
for i = 1:p
    z = [z, phif(:, :, i+1)];
end
for i = 1:q
    z = [z, thf(:, :, i+1)];
end
z = [z, Sigmar];
mprint(z, in);
disp('press any key to continue')
pause

%residual diagnostics
e = result.e; %white noise residuals
ff = result.ff; %vector of nonlinear functions
nbeta = 0; %length of regression vector
ndrs = length(yd);
Ss = e' * e;
Ff = ff' * ff; %residual sum of squares
conp = result.sigma2c;
sconp = sqrt(conp);
lagl = 26;
infr = rescomp(e, lagl, length(xvf), Ss, conp, sconp, Ff, ndrs, nbeta);

ic = 1;
lag = 26;
disp(' ')
disp('******** Residuals:     ********');
str = mautcov(e, lag, ic);
disp('More estimation and diagnostic details are in file: ');
disp(tname)
disp('in the subdirectory "results"')
disp('press any key to continue')
pause

%plot residual diagnostics
plotres([], [], [], [], [], 1.96, 'residuals', 1, 0, [], 0, [], infr, 1, 1);
close all

%print residual diagnostics
%file for output
fname = ['results', filesep, tname, '.txt'];
fid = fopen(fname, 'w');
% fid=1;
printres(fid, infr);
%close external file
if fid ~= 1
    fclose(fid);
end


%compute forecasts of the original series
if mpr > 0
    npr = mpr;
    %set up system matrices for the estimated ARIMA model
    %Differencing polynomial
    diff(:, :, 1) = 1.;
    if (d > 0)
        diff1(:, :, 1) = 1.;
        diff1(:, :, 2) = -1.;
        for i = 1:d
            diff = pmatmul(diff, diff1);
        end
    end
    %Regular AR polynomial for the original series: the product of diff and
    %phif
    [phio, ierror] = pmatmul(phif, diff);
    %Regular MA polynomial as before
    tho = thf;
    %Seasonal MA polynomials as before
    Phio = Phif;
    Tho = Thf;
    %Note that the residual covariance matrix is divided by the concentrated
    %parameter (result.sigma2c).
    Sigmaf = Lf * Lf';
    %set up state space model
    [stro, ferror] = suvarmapqPQ(phio, tho, Phio, Tho, Sigmaf, freq);
    %ARIMA model in state space form
    Z = stro.Z;
    G = stro.G;
    T = stro.T;
    H = stro.H;
    [ndelta, junk] = size(T);
    X = Y;
    W = [];
    %initial conditions for the Kalman filter. They are obtained automatically
    %by the following function.
    ndelta = d; %number of unit roots in the model
    [ins, i, ferror] = incossm(T, H, ndelta);
    if isempty(Y)
        chb = 0; %there are no regression effects, so do not compute hb and Mb in
        %scakfle2
    else
        chb = 1;
    end
    
    %run Kalman filter
    [e, f, hb, Mb, A, P, qyy, R] = scakfle2(y, X, Z, G, W, T, H, ins, i, chb);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X = [];
    chb = 0;
    [e, f, hb, Mb, A, P, qyy, R] = scakfle2(y, X, Z, G, W, T, H, ins, i, chb);
    [KKP, PT, hd, Md, initf, recrs, recr, srecr] = scakff(y, X, Z, G, W, T, H, ins, i);
    return
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %hb is the vector of regression estimates and Mb is the matrix of standard
    %errors. A is the estimated state vector, x_{t|t-1}, obtained with the
    %Kalman filter at the end of the sample and P is the matrix of standard
    %errors.
    
    %forecasts
    m = 1; %number of series
    if isempty(Y)
        Xp = Y;
    else
        Xp = YY(n+1:n+npr, :);
    end
    Wp = [];
    cw = 1.96;
    [pry, mypr, alpr, malpr] = ssmpred(npr, m, A, P, Xp, Z, G, Wp, T, H, hb, Mb);
    spry = zeros(m, npr);
    sconp = sqrt(result.sigma2c);
    for i = 1:npr
        spry(:, i) = sqrt(diag(mypr(:, :, i))) * sconp;
    end
    opry = pry;
    ospry = spry;
    %plot forecasts
    out.pry = pry;
    out.spry = spry;
    out.opry = opry;
    out.ospry = ospry;
    out.y = y;
    out.yor = y;
    out.ny = length(y);
    out.npr = npr;
    out.cw = cw;
    out.tname = tname;
    out.lam = lam;
    out.s = freq;
    pfctsusm(out);
end
