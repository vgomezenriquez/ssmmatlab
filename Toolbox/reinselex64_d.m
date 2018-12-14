%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example 6.4 of Reinsel (1997), pp. 206-209
%
% Series are: 1) US housing starts and 2) US housing sold.
% The period is: January 1965 through December 1974.
% The number of observations is ny=120.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;

y = load(fullfile('data', 'housingt.dat'));
[ny, my] = size(y);

lag = 36;
cw = 1.96;
freq = 12;
dr = 0;
tname = {'US housing starts', 'US housing sold'};
for i = 1:2
    for ds = 0:1
        c0 = sacspacdif(y(:, i), tname(i), dr, ds, freq, lag, cw);
        disp('press any key to continue')
        pause
    end
end
close all

% %identify a VAR model for the original series
% maxlag=6; minlag=1; prt=1;
% lagsopt = varident(y,maxlag,minlag,prt);
% disp('estimated VAR lag length: ')
% disp(lagsopt)


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
%estimate ''differencing polynomial'' D and obtain differenced series yd
%note that betaor is parametrized in DA
[D, DA, yd, ferror] = mdfestim1r(y, x, prt, nr);
seas = 1;

disp(' ')
disp('identify a VAR model for the differenced series')
disp(' ')
disp('press any key to continue')
pause
%identify a VAR model for the differenced series
maxlag = 3;
minlag = 1;
prt = 1;
lagsopt = varident(yd, maxlag, minlag, prt);
pause

disp(' ')
disp('estimation of a VAR(1) model for the differenced series')
disp('')
disp('press any key to continue')
pause

%estimation of a VAR(1) model for the differenced series
hr3 = 0;
finv2 = 1;
mstainv = 1;
pr = 1;
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
% Example 6.4
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
disp('estimation of the VAR(1) model in error correction form')
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
mprintar(Lambdaf, in, tit);
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
%
% End of Example 6.4
% Estimation of the model parameterized in error correction form
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example 6.4
% Estimation of the model parameterized in terms of the model for the
% differenced series
%

%set up varma model for the differenced series
clear str
[str, ferror] = suvarmapqPQ(phi, th, Phi, Th, Sigma, seas);

%add unit root information to the model structure
[str, ferror] = aurirvarmapqPQ(str, nr, DA);

disp(' ')
disp('estimation of the model in terms of the model for ')
disp('the differenced series. Note that the overall model')
disp('has changed. It is now VAR(2).                     ')
disp('press any key to continue')
pause

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
disp(' ');
disp('***** Estimated Model  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'Lambda';
strt = 1;
mprintar(Lambda(:, :, 2), in, tit, strt);
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
%
% end of Example 6.4
% Estimation of the model parameterized in terms of the model for the
% differenced series
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
