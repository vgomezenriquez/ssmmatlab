%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example 6.1 of Tsay (2014), pp. 335-341
%
% Monthly housing data of the United States from
% January 1963 to July 2012. The two series employed are
% 1. z1t: Logarithm of new homes sold in thousands of units (new residential sales)
% 2. z2t: Logarithm of the total new privately owned housing units started in
%             thousands of units (new residential construction)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

da = load(fullfile('data', 'm-hsoldhst6312.dat'));
zt = log(da(:, 3:4));
tdx = da(:, 1) + da(:, 2) / 12;
[nx, mx] = size(zt);
freq = 12;

subplot(2, 1, 1)
plot(tdx, zt(:, 1))
xlabel('time');
ylabel('Hsold');
axis('tight');
subplot(2, 1, 2)
plot(tdx, zt(:, 2))
xlabel('time');
ylabel('Hstart');
axis('tight');
disp('press any key to continue')
pause
close all


lag = 36;
cw = 1.96;
dr = 0;
tname = {'US housing sold', 'US housing starts'};
for i = 1:2
    for ds = 0:1
        c0 = sacspacdif(zt(:, i), tname(i), dr, ds, freq, lag, cw);
        pause
    end
end
close all


yd = diferm(zt, freq); %seasonal differencing
yd = diferm(yd, 1); %regular differencing

%compute autocovariance and autocorrelation matrices
lag = 24;
ic = 1;
nr = 0;
disp(' ')
disp('******** Sample cross correlation matrics:     ********');
stre = mautcov(yd, lag, ic, nr);
disp('Correlation matrix at lag 0:')
disp(stre.r0)


disp(' ')
disp('estimate VARMA(0,3)(0,1) model using the Hannan-Rissanen method ')
disp('and eliminate some insignificant parameters')
disp('press any key to continue')
pause


%estimate a VARMAX(0,3,0)(0,1,0) model by the Hannan-Rissanen method.
x = [];
hr3 = 0;
finv2 = 1;
mstainv = 1;
nsig = [1, 1];
tsig = [.8, .8];
qr = 3;
qs = 1;
[strv, ferror] = estvarmaxpqrPQR(yd, x, freq, [0, qr, 0], [0, qs, 0], hr3, finv2, mstainv, nsig, tsig);

disp(' ');
disp('***** Estimated VARMA(0,3)(0,1) Model using the HR method  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'th';
strt = 1;
mprintar(strv.thetas3(:, :, 2:qr+1), in, tit, strt);
disp(' ');
% in.cnames = char('Th(1):',' ');
% mprint(strv.thetas3(:,:,freq+qs),in);
tit = 'Th';
strt = 1;
mprintar(strv.thetas3(:, :, freq+qs), in, tit, strt);
disp('press any key to continue')
pause


%estimate using exact ML
%setup model
Phi = eye(2);
Th(:, :, 1) = eye(2);
Th(:, :, qs+1) = strv.thetas3(:, :, freq+qs);
th = strv.thetas3(:, :, 1:qr+1);
phi = eye(2);
Sigma = strv.sigmar3;

[str, ferror] = suvarmapqPQ(phi, th, Phi, Th, Sigma, freq);

%eliminate insignificant parameters
for k = 2:qr + 1
    for i = 1:mx
        for j = 1:mx
            if th(i, j, k) == 0
                str.thn(i, j, k) = 0;
            end
        end
    end
end
for i = 1:mx
    for j = 1:mx
        if Th(i, j, 2) == 0
            str.Thn(i, j, 2) = 0;
        end
    end
end
[str, ferror] = fixvarmapqPQ(str);


Y = [];

disp(' ')
disp('estimation using the exact method')
disp('press any key to continue')
pause

%estimate model using the exact method
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

%create estimated model
[phif, thf, Phif, Thf, Lf, ferror] = pr2varmapqPQ(xvf, xf, str);
Sigmar = result.Sigmar;
%t-values
tvf = result.tv;


[phitvf, thtvf, Phitvf, Thtvf, Ltvf, ferror] = pr2varmapqPQ(tvf, xf, str);

disp(' ');
disp('***** Estimated Model  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'th';
strt = 1;
mprintar(thf(:, :, 2:qr+1), in, tit, strt);
disp(' ')
% in.cnames = char('Th(1):',' ');
% mprint(Thf(:,:,2),in);
tit = 'Th';
strt = 1;
mprintar(Thf(:, :, 2), in, tit, strt);
tit = 'Sigma:';
mprintar(Sigmar, in, tit);
disp(' ')
disp('press any key to continue')
pause

disp(' ');
disp('***** t-values  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'tv-th';
strt = 1;
mprintar(thtvf(:, :, 2:qr+1), in, tit, strt);
disp(' ')
% in.cnames = char('tv-Th(1):',' ');
% mprint(Thtvf(:,:,2),in);
tit = 'tv-Th';
strt = 1;
mprintar(Thtvf(:, :, 2), in, tit, strt);
% in.cnames = char('tv-L:',' ');
% mprint(Ltvf,in);
tit = 'tv-L';
mprintar(Ltvf, in, tit);
disp('press any key to continue')
pause

disp(' ')
disp('******** Computation of Recursive Residuals   ********');
%compute recursive residuals
[strf, ferror] = suvarmapqPQ(phif, thf, Phif, Thf, Sigmar, freq);
%set up regression matrices
X = Y;
W = [];
%set up system matrices
T = strf.T;
Z = strf.Z;
G = strf.G;
H = strf.H;
%set up initial conditions
ndelta = 0; %number of unit roots
[ins, i, ferror] = incossm(T, H, ndelta);

[Xt, Pt, g, M, initf, recrs, recr] = scakff(yd, X, Z, G, W, T, H, ins, i);
%plot recursive residuals
plot(recr(:, 1)), legend('Hsold'), pause
plot(recr(:, 2)), legend('Hstart'), pause
close all
%compute autocovariance and autocorrelation matrices of rec. residuals
lag = 24;
ic = 1;
nr = length(result.xvf);
disp(' ')
disp('******** Recursive Residuals:     ********');
str = mautcov(recr, lag, ic, nr);
disp('Correlation matrix at lag 0:')
disp(str.r0)
disp('Q statistics:')
disp(str.qstat)

disp('p-values of Q statistics:')
disp(str.pval)
[m, n] = size(str.pval);
t = 1:m;
plot(t, str.pval, t, 0.05*ones(1, m))
legend('p-values of Q statistics:')
pause
close all

%obtain 24 forecasts for the differenced series
mpr = 24;
%hb, Mb, A and P are in structure result. Here, hb is the vector of
%regression estimates and Mb is the matrix of mse. A is the
%estimated state vector, x_{t|t-1}, obtained with the Kalman filter at the
%end of the sample and P is the matrix of mse.
hb = result.h;
Mb = result.H;
A = result.A;
P = result.P;

npr = mpr; %number of forecasts
Xp = Y;
Wp = [];
m = 2; %number of series
[pry, mypr, alpr, malpr] = ssmpred(npr, m, A, P, Xp, Z, G, Wp, T, H, hb, Mb);
conp = result.sigma2c;
spry = zeros(m, npr);
sconp = sqrt(result.sigma2c);
for i = 1:npr
    spry(:, i) = sqrt(diag(mypr(:, :, i))) * sconp;
end
opry = pry;
ospry = spry;
%plot forecasts
tname = 'var1';
out.pry = pry(1, :);
out.spry = spry(1, :);
out.opry = opry(1, :);
out.ospry = ospry(1, :);
out.y = yd(:, 1);
out.yor = yd(:, 1);
out.ny = length(yd(:, 1));
out.npr = npr;
out.cw = cw;
out.tname = tname;
lam = 1; %lam=0, logs are taken; =1, no logs are taken
out.lam = lam;
out.s = freq;
pfctsusm(out);
tname = 'var2';
out.pry = pry(2, :);
out.spry = spry(2, :);
out.opry = opry(2, :);
out.ospry = ospry(2, :);
out.y = yd(:, 2);
out.yor = yd(:, 2);
out.ny = length(yd(:, 2));
out.npr = npr;
out.cw = cw;
out.tname = tname;
lam = 1; %lam=0, logs are taken; =1, no logs are taken
out.lam = lam;
out.s = freq;
pfctsusm(out);

%compute 24 forecasts of the original series
%set up system matrices for the estimated VARMA model, including the
%differencing matrix polynomial.
%Differencing polynomial
phifo(:, :, 1) = eye(2);
phifo(:, :, 2) = -eye(2);
Phifo(:, :, 1) = eye(2);
Phifo(:, :, 2) = -eye(2);
%MA polynomial
thfo = thf;
Thfo = Thf;
[strfo, ferror] = suvarmapqPQ(phifo, thf, Phifo, Thf, Sigmar, freq);
%ARIMA model in state space form
Z = strfo.Z;
G = strfo.G;
T = strfo.T;
H = strfo.H;
[ndelta, junk] = size(T);
X = [];
W = [];
%initial conditions for the Kalman filter
[ins, ii, ferror] = incossm(T, H, ndelta);
chb = 0; %there are no regression effects, so do not compute hb and Mb in
%scakfle2

%run Kalman filter
[e, f, hb, Mb, A, P, qyy, R] = scakfle2(zt, X, Z, G, W, T, H, ins, ii, chb);
%hb is the vector of regression estimates and Mb is the matrix of standard
%errors. A is the estimated state vector, x_{t|t-1}, obtained with the
%Kalman filter at the end of the sample and P is the matrix of standard
%errors.

%forecasts
[pry, mypr, alpr, malpr] = ssmpred(npr, m, A, P, Xp, Z, G, Wp, T, H, hb, Mb);
spry = zeros(m, npr);
%result.sigma2c is the (1,1) parameter in the covariance matrix of the
%innovations, that is always concentrated out of the likelihood
sconp = sqrt(result.sigma2c);
for i = 1:npr
    spry(:, i) = sqrt(diag(mypr(:, :, i))) * sconp;
end
opry = pry;
ospry = spry;
%plot forecasts
tname = 'var1';
out.pry = pry(1, :);
out.spry = spry(1, :);
out.opry = opry(1, :);
out.ospry = ospry(1, :);
out.y = zt(:, 1);
out.yor = zt(:, 1);
out.ny = length(zt(:, 1));
out.npr = npr;
out.cw = cw;
out.tname = tname;
lam = 1; %lam=0, logs are taken; =1, no logs are taken
out.lam = lam;
out.s = freq;
pfctsusm(out);
tname = 'var2';
out.pry = pry(2, :);
out.spry = spry(2, :);
out.opry = opry(2, :);
out.ospry = ospry(2, :);
out.y = zt(:, 2);
out.yor = zt(:, 2);
out.ny = length(zt(:, 2));
out.npr = npr;
out.cw = cw;
out.tname = tname;
lam = 1; %lam=0, logs are taken; =1, no logs are taken
out.lam = lam;
out.s = freq;
pfctsusm(out);
