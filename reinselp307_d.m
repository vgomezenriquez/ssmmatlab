%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Data given in Reinsel (1997), p. 307
%
% Three US monthly interest rates from 1960 to 1979 (Stock and Watson,
% 1988; Reinsel, 1997, p. 307). The length is n=240. Yap and Reinsel (1995)
% estimate a VARMA(1,1) with cointegration rank equal to two. Stock and
% Watson (1988) also analyzed this series.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;

yy = load(fullfile('data', 'monthlyusir.dat'));
y = yy(:, 3:end);
yor = y;
y = log(y);

[ny, my] = size(y);

% lag=36; cw=1.96; freq=1; ds=0;
% tname={'Series 1','Series 2','Series 3'};
% for i=1:my
%  for dr=0:1
%   c0=sacspacdif(y(:,i),tname(i),dr,ds,freq,lag,cw);
%   pause
%  end
% end
% close all

disp(' ')
disp('identify a VAR model for the original series')
disp(' ')
disp('press any key to continue')
pause
%identify a VAR model for the original series
maxlag = 6;
minlag = 1;
prt = 1;
lagsopt = varident(y, maxlag, minlag, prt);
disp('press any key to continue')
pause

% disp(' ')
% disp('identify a VARMA(p,q) model for the original series')
% disp(' ')
% disp('press any key to continue')
% pause
%
% %identify a VARMA(p,q) model for the series
% maxlag=6; minlag=0; prt=0; x=[]; seas=1;
% [lagsopt,ferror] = lratiopqr(y,x,seas,maxlag,minlag,prt);
% disp(' ')
% disp('Estimated orders in VARMAX(p,q,r):  ')
% disp(lagsopt)
% disp('press any key to continue')
% pause


% %Matlab-Econ function for Johansen cointegration test
% lags=3;
% [h,pValue,stat,cValue,mles] = jcitest(y,'lags',lags);pause

[D, nr, yd, DA, ferror] = mcrcregr(y);
disp(' ')
disp('number of unit roots according to the crc criterion:')
disp(nr)
disp('press any key to continue')
pause


%number of unit roots in the model
nr = 1;
prt = 0;
x = [];
%estimate ''differencing polynomial'' D and obtain differenced series yd
%note that betaor is parametrized in DA
[D, DA, yd, ferror] = mdfestim1r(y, x, prt, nr);
seas = 1;


disp(' ')
disp('estimation of a VARMA(1,1) model for the original series')
disp('')
disp('press any key to continue')
pause

%estimation of the model for the original series
% model for the original series: VARMA(1,1)
hr3 = 0;
finv2 = 1;
mstainv = 1;
pr = 1;
qr = 1; %VARMA(1,1)
[str, ferror] = estvarmaxpqrPQR(y, x, seas, [pr, qr, 0], [0, 0, 0], hr3, finv2, mstainv);

phi(:, :, 1) = str.phis3(:, :, 1);
for i = 1:pr
    phi(:, :, i+1) = str.phis3(:, :, i+1);
end
Phi(:, :, 1) = phi(:, :, 1);
th(:, :, 1) = str.thetas3(:, :, 1);
for i = 1:qr
    th(:, :, i+1) = str.thetas3(:, :, i+1);
end
Th(:, :, 1) = phi(:, :, 1);
Sigma = str.sigmar3;
clear str


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example of Yap and Reinsel (1995)
% Estimation of the model parameterized in in error correction form
%

%pass AR part to error correction form
phit = pmatmul(phi, Phi);
[Pi, Lambdap, alpha, betap, ferror] = mid2mecf(phit, D, DA);
disp(' ');
disp('***** Estimated Model  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
C = -eye(my) - phi(:, :, 2);
tit = 'C';
mprintar(C, in, tit);
tit = 'Theta';
mprintar(th(:, :, 2), in, tit);
disp('press any key to continue')
pause


%set up varma model in error correction form. Parameters are defined in
%terms of the error correction form.
%re-specify Lambda so that the model coincides with that of Yap and Reinsel
%(1995)
Lambda = Lambdap(:, :, 1);
[str, ferror] = suvarmapqPQe(Lambda, alpha, betap, th, Th, Sigma, seas);
return


disp(' ')
disp('estimation of the VARMA(1,1) model in error correction form')
disp('press any key to continue')
pause

%estimate the model in error correction form
Y = [];
constant = 1;
result = varmapqPQestime(y, str, Y, constant);

disp(' ');
disp('******************** Results from estimation ********************');
disp(' ');
mprintr(result)
disp('press any key to continue')
pause


%create estimated model
xvrf = result.xvf;
xrf = result.xf;
[ydf, xvv, xff, DAf, Dr, Ds, ferror] = pr2varmapqPQd(y, xvrf, xrf, str);
[Lambdaf, alphaf, betapf, thf, Thf, Lf, ferror] = pr2ecf(xvv, xff, DA, str);
[phif, D, DA, ferror] = mecf2mid(Lambdaf, alphaf, betapf);
disp(' ');
disp('***** Estimated Model  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'Lambda';
mprintar(Lambdaf(:, :, 1), in, tit);
tit = 'Theta';
strt = 1;
mprintar(thf(:, :, 2), in, tit, strt);
%reparameterize betap and alpha
[nb, mb] = size(betapf);
alphaf = alphaf * betapf(:, 1:nb);
betapf = betapf(:, 1:nb) \ betapf;
disp('estimated beta transposed')
disp(betapf)
disp('estimated alpha')
disp(alphaf)
%constant is phif(1)*mu
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'Sigma:';
mprintar(result.Sigmar, in, tit);
tit = 'Constant:';
mprintar((sum(phif, 3) * result.h)', in, tit);
disp('press any key to continue')
pause

%
% End of Example of Yap and Reinsel (1995)
% Estimation of the model parameterized in error correction form
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%forecasting
%form overall AR matrix polynomial
phit = pmatmul(phif, D);
%last matrix in phit is insignificant
phit = phit(:, :, 1:end-1);
%create structure and put model into state space form
Sigma = result.Sigmar;
[str, ferror] = suvarmapqPQ(phit, thf, Phi, Th, Sigma, seas);
%model matrices
T = str.T;
Z = str.Z;
G = str.G;
H = str.H;
W = str.W;
%initial conditions
%we have to pass the number of unit roots in the following function to
%obtain the initial conditions.
[ins, ii, ferror] = incossm(T, H, nr);

%number of forecasts
npr = 24;

%generate regression variables for the original series. Since
%D(B)Y_t = c + U_t,
%where D(z) is the differencing matrix polynomial, the variable is
%X_t = D^{-1}(B)I, or D(B)*X_t = I.
%
nypnpr = ny + npr;
nx = nypnpr * my;
X = zeros(nx, my);
for i = 2:nypnpr
    ipm1 = (i - 2) * my + 1:(i - 1) * my;
    ip = (i - 1) * my + 1:i * my;
    X(ip, :) = eye(my) - D(:, :, 2) * X(ipm1, :);
end
%run Kalman filter
nx1 = ny * my;
X1 = X(1:nx1, :);
chb = 1;
[e, f, hb, Mb, A, P, qyy, R] = scakfle2(y, X1, Z, G, W, T, H, ins, ii, chb);
%hb is the vector of regression estimates and Mb is the matrix of standard
%errors. A is the estimated state vector, x_{t|t-1}, obtained with the
%Kalman filter at the end of the sample and P is the matrix of standard
%errors.

%forecasts in logs
%regression variables for the forecasts
Xp = X(nx1+1:nx, :);
[pry, mypr, alpr, malpr] = ssmpred(npr, my, A, P, Xp, Z, G, W, T, H, hb, Mb);
spry = zeros(my, npr);
sconp = sqrt(result.sigma2c);
for i = 1:npr
    spry(:, i) = sqrt(diag(mypr(:, :, i))) * sconp;
end
cw = 1.96;
lam = 1;
opry = pry;
ospry = spry;
%plot forecasts
for i = 1:my
    tname = ['var', num2str(i), ' in logs'];
    out.pry = pry(i, :);
    out.spry = spry(i, :);
    out.opry = opry(i, :);
    out.ospry = ospry(i, :);
    out.y = y(:, i);
    out.yor = y(:, i);
    out.ny = length(y(:, i));
    out.npr = npr;
    out.cw = cw;
    out.tname = tname;
    out.lam = lam;
    out.s = seas;
    pfctsusm(out);
    clear out
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example of Yap and Reinsel (1995)
% Estimation of the VARMA(1,1) model parameterized in terms of the model
% for the differenced series.
%

disp(' ')
disp('estimation of the VARMA(1,1) model for the differenced series.')
disp('Note that the overall model has changed. It is now VARMA(2,1).')
disp('press any key to continue')
pause


%set up varma model for the differenced series. Parameters are defined in
%terms of the model for the differenced series, contained in phi, th, Phi,
%Th and Sigma, and the parameters of the differencing matrix polynomial,
%contained in DA.
clear str
[str, ferror] = suvarmapqPQ(phi, th, Phi, Th, Sigma, seas);

%add unit root information to the model structure
[str, ferror] = aurirvarmapqPQ(str, nr, DA);

%estimate model with cointegration relations imposed
Y = [];
constant = 1;
result = varmapqPQestimd(y, str, Y, constant);

disp(' ');
disp('******************** Results from estimation ********************');
disp(' ');
mprintr(result)
disp('press any key to continue')
pause

%create estimated model
xvrf = result.xvf;
xrf = result.xf;
[ydf, xvv, xff, DAf, Dr, Ds, ferror] = pr2varmapqPQd(y, xvrf, xrf, str);
[phif, thf, Phif, Thf, Lf, ferror] = pr2varmapqPQ(xvv, xff, str);


%error correction matrices
[Pi, Lambda, alpha, betafp, ferror] = mid2mecf(phif, Dr, DAf);
%note that the Lambda matrix polynomial has degree one in this model
%whereas it had degree zero in the previous model.
disp(' ');
disp('***** Estimated Model  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'Lambda';
strt = 1;
mprintar(Lambda(:, :, 2), in, tit, strt);
tit = 'Theta';
strt = 1;
mprintar(thf(:, :, 2), in, tit, strt);
%reparameterize betap and alpha
[nb, mb] = size(betafp);
alpha = alpha * betafp(:, 1:nb);
betafp = betafp(:, 1:nb) \ betafp;
disp('estimated beta transposed')
disp(betafp)
disp('estimated alpha')
disp(alpha)
%constant is phif(1)*mu
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'Sigma:';
mprintar(result.Sigmar, in, tit);
tit = 'Constant:';
mprintar((sum(phif, 3) * result.h)', in, tit);
disp('press any key to continue')
pause

disp(' ')
disp('Insignificant parameters are fixed to zero ')
disp('and the model is reestimated')
disp('press any key to continue')
pause

%fix insignificant parameters in phif to zero using the field phin in a
%new structure str and estimate model again
phif(3, 2, 2) = 0.0;
phif(2, 3, 2) = 0.0;
thf(2, 1, 2) = 0.0;
thf(3, 1, 2) = 0.0;
% thf(2,2,2)=0.0; thf(3,2,2)=0.0;
clear str
[str, ferror] = suvarmapqPQ(phif, thf, Phif, Thf, result.Sigmar, seas);
%add unit root information to the model structure
[str, ferror] = aurirvarmapqPQ(str, nr, DAf);
str.phin(3, 2, 2) = 0.0;
str.phin(2, 3, 2) = 0.0;
str.thn(2, 1, 2) = 0.0;
str.thn(3, 1, 2) = 0.0;
[str, ferror] = fixvarmapqPQ(str);

%estimate model again
resultf = varmapqPQestimd(y, str, Y, constant);

disp(' ');
disp('******************** Results from estimation ********************');
disp(' ');
mprintr(resultf)
disp('press any key to continue')
pause


%create new estimated model
xvrff = resultf.xvf;
xrff = resultf.xf;
[ydff, xvvf, xfff, DAff, Drf, Dsf, ferror] = pr2varmapqPQd(y, xvrff, xrff, str);
[phiff, thff, Phiff, Thff, Lff, ferror] = pr2varmapqPQ(xvvf, xfff, str);

disp('estimated residual covariance matrix')
disp(resultf.Sigmar)


%error correction matrices
[Pif, Lambdaf, alphaf, betafpf, ferror] = mid2mecf(phiff, Drf, DAff);

disp(' ');
disp('***** Estimated Model  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'Lambda';
strt = 1;
mprintar(Lambdaf(:, :, 2), in, tit, strt);
tit = 'Theta';
strt = 1;
mprintar(thff(:, :, 2), in, tit, strt);
%reparameterize betap and alpha
[nb, mb] = size(betafpf);
alphaf = alphaf * betafpf(:, 1:nb);
betafpf = betafpf(:, 1:nb) \ betafpf;
disp('estimated beta transposed')
disp(betafpf)
disp('estimated alpha')
disp(alphaf)
%constant is phif(1)*mu
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'Sigma:';
mprintar(resultf.Sigmar, in, tit);
tit = 'Constant:';
mprintar((sum(phiff, 3) * result.h)', in, tit);
%
% end of Example of Yap and Reinsel (1995)
% Estimation of the model parameterized in terms of the model for the
% differenced series
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
