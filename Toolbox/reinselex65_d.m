%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example 6.5 of Reinsel (1997), pp. 213-215
%
% Series are: U.S. quarterly interest rate on 1) AAA corporate bonds and on
% 2) commercial paper.
% The number of observations is ny=72. A VAR(3) model with cointegration
% rank equal to one is estimated.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;

yy = load(fullfile('data', 'aaa-cp.dat'));
y = yy(:, 3:4);

[ny, my] = size(y);

lag = 20;
cw = 1.96;
freq = 1;
ds = 0;
tname = {'interest rate on AAA corporate bonds', ...
    'interest rate on commercial paper'};
for i = 1:2
    for dr = 0:1
        c0 = sacspacdif(y(:, i), tname(i), dr, ds, freq, lag, cw);
        pause
    end
end
close all


% %Matlab-Econ function for Johansen cointegration test
% lags=lagsopt;
% [h,pValue,stat,cValue,mles] = jcitest(y,'lags',lags);

[D, nr, yd, DA, ferror] = mcrcregr(y);
disp('number of unit roots according to the crc criterion:')
disp(nr)
disp('press any key to continue')
pause

%number of unit roots in the model
nr = 1;
prt = 0;
x = [];
disp(' ')
disp('estimate ''differencing polynomial'' and obtain differenced series')
disp(' ')
disp('press any key to continue')
pause
%estimate ''differencing matrix polynomial'' D and obtain differenced
%series yd.
%note that betaor is parametrized in DA
[D, DA, yd, ferror] = mdfestim1r(y, x, prt, nr);

%identify a VAR model for the differenced series
maxlag = 6;
minlag = 1;
prt = 1;
lagsopt = varident(yd, maxlag, minlag, prt);
disp('estimated VAR lag length for the differenced series: ')
disp(lagsopt)
disp('press any key to continue')
pause

disp(' ')
disp('estimation of a VAR(3) model for the differenced series')
disp('')
disp('press any key to continue')
pause


%estimation of a VAR(3) model for the differenced series
seas = 1;
hr3 = 0;
finv2 = 1;
mstainv = 1;
pr = 3;
[str, ferror] = estvarmaxpqrPQR(yd, x, seas, [pr, 0, 0], [0, 0, 0], hr3, finv2, mstainv);

phi(:, :, 1) = str.phis3(:, :, 1);
for i = 2:pr + 1
    phi(:, :, i) = str.phis3(:, :, i);
end
Phi(:, :, 1) = phi(:, :, 1);
th(:, :, 1) = phi(:, :, 1);
Th(:, :, 1) = phi(:, :, 1);
Sigma = str.sigmar3;
clear str


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example 6.5
% Estimation of the model parameterized in in error correction form
%

%pass AR part to error correction form
phit = pmatmul(phi, Phi);
[Pi, Lambda, alpha, betap, ferror] = mid2mecf(phit, D, DA);
%re-specify Lambda so that the model coincides with that of Reinsel (1997)
Lambda = Lambda(:, :, 1:end-1);

%set up varma model in error correction form. Parameters are defined in
%terms of the error correction form.
[str, ferror] = suvarmapqPQe(Lambda, alpha, betap, th, Th, Sigma, seas);

disp(' ')
disp('estimation of the VAR(3) model in error correction form')
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
disp(' ');
disp('***** Estimated Model  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'Lambda';
strt = 1;
mprintar(Lambdaf(:, :, 2:3), in, tit, strt);
disp(' ')
%reparameterize betap and alpha
[nb, mb] = size(betapf);
alphaf = alphaf * betapf(:, 1:nb);
betapf = betapf(:, 1:nb) \ betapf;
disp('estimated beta transposed')
disp(betapf)
disp('estimated alpha')
disp(alphaf)
disp('estimated residual covariance matrix')
disp(result.Sigmar)

disp(' ')
disp('Insignificant parameters are fixed to zero ')
disp('and the model is reestimated')
disp('press any key to continue')
pause


%fix insignificant parameters Lambda(1,2,2), Lambda(1,2,3) and
%Lambda(1,3,3) to zero using the field Lambdan
clear str
Lambdaf(1, 2, 2) = 0.0;
Lambdaf(1, 2, 3) = 0.0;
Lambdaf(2, 2, 3) = 0.0;
[str, ferror] = suvarmapqPQe(Lambdaf, alphaf, betapf, thf, Thf, result.Sigmar, seas);
str.Lambdan(1, 2, 2) = 0.0;
str.Lambdan(1, 2, 3) = 0.0;
str.Lambdan(2, 2, 3) = 0.0;
%if we do not fix this parameter to its actual value, it does not converge
str.Lhn(2) = str.Lh(2);
[str, ferror] = fixvarmapqPQe(str);
%estimate model again
resultf = varmapqPQestime(y, str, Y, constant);

disp(' ');
disp('******************** Results from estimation ********************');
disp(' ');
mprintr(resultf)
disp('press any key to continue')
pause


%create estimated model
xvrf = resultf.xvf;
xrf = resultf.xf;
[ydf, xvv, xff, DAff, Dr, Ds, ferror] = pr2varmapqPQd(y, xvrf, xrf, str);
[Lambdaff, alphaff, betapff, thff, Thff, Lff, ferror] = pr2ecf(xvv, xff, DAff, str);
disp(' ');
disp('***** Estimated Model  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'Lambda';
strt = 1;
mprintar(Lambdaff(:, :, 2:3), in, tit, strt);
disp(' ')
%reparameterize betap and alpha
[nb, mb] = size(betapff);
alphaff = alphaff * betapff(:, 1:nb);
betapff = betapff(:, 1:nb) \ betapff;
disp('estimated beta transposed')
disp(betapff)
disp('estimated alpha')
disp(alphaff)
disp('estimated residual covariance matrix')
disp(resultf.Sigmar)
%
% End of Example 6.5
% Estimation of the model parameterized in error correction form
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example 6.5
% Estimation of the model parameterized in terms of the model for the
% differenced series
%

disp(' ')
disp('estimation of the model in terms of the model for ')
disp('the differenced series. Note that the overall model')
disp('has changed. It is now VAR(4).                     ')
disp('press any key to continue')
pause


%set up varma model for the differenced series. Parameters are defined in
%terms of the model for the differenced series, contained in phi, th, Phi,
%Th and Sigma, and the parameters of the differencing matrix polynomial,
%contained in DA.
[str, ferror] = suvarmapqPQ(phi, th, Phi, Th, Sigma, seas);


%add unit root information to the model structure
[str, ferror] = aurirvarmapqPQ(str, nr, DA);

Y = [];
constant = 0;
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
disp(' ');
disp('***** Estimated Model  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'Phi';
strt = 1;
mprintar(Lambda(:, :, 2:4), in, tit, strt);
disp(' ')
%reparameterize betap and alpha
[nb, mb] = size(betafp);
alpha = alpha * betafp(:, 1:nb);
betafp = betafp(:, 1:nb) \ betafp;
disp('estimated hat B')
disp(betafp)
disp('estimated hat A')
disp(alpha)
disp('estimated residual covariance matrix')
disp(result.Sigmar)

disp(' ')
disp('Insignificant parameters are fixed to zero ')
disp('and the model is reestimated')
disp('press any key to continue')
pause

%fix insignificant parameter phif(2,2,3) to zero using the field phin in a
%new structure str and estimate model again
phif(2, 2, 3) = 0.0;
clear str
[str, ferror] = suvarmapqPQ(phif, thf, Phif, Thf, result.Sigmar, seas);
%add unit root information to the model structure
[str, ferror] = aurirvarmapqPQ(str, nr, DAf);
str.phin(2, 2, 3) = 0.0;
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


%error correction matrices
[Pif, Lambdaf, alphaf, betafpf, ferror] = mid2mecf(phiff, Drf, DAff);
disp(' ');
disp('***** Estimated Model  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'Phi';
strt = 1;
mprintar(Lambdaf(:, :, 2:4), in, tit, strt);
disp(' ')
%reparameterize betap and alpha
[nb, mb] = size(betafpf);
alphaf = alphaf * betafpf(:, 1:nb);
betafpf = betafpf(:, 1:nb) \ betafpf;
disp('estimated hat B')
disp(betafpf)
disp('estimated hat A')
disp(alphaf)
disp('estimated residual covariance matrix')
disp(resultf.Sigmar)
%
% end of Example 6.5
% Estimation of the model parameterized in terms of the model for the
% differenced series
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
