%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example 6.6 of Reinsel (1997), pp. 221-224
%
% Series are: 1) US housing starts and 2) US housing sold.
% The period is: January 1965 through December 1974.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;

yt = load(fullfile('data', 'housing.dat'));
y = yt(:, 3:4);

lag = 36;
cw = 1.96;
freq = 12;
dr = 0;
tname = {'US housing starts', 'US housing sold'};
for i = 1:2
    for ds = 0:1
        c0 = sacspacdif(y(:, i), tname(i), dr, ds, freq, lag, cw);
        pause
    end
end
close all


yd = diferm(y, freq);

disp(' ')
disp('estimate VARMA(1,0)(0,1) model using the Hannan-Rissanen method: ')
disp('press any key to continue')
pause


%estimate a VARMA(1,0,0)(0,1,0) model by the Hannan-Rissanen method.
x = [];
hr3 = 0;
finv2 = 1;
[strv, ferror] = estvarmaxpqrPQR(yd, x, freq, [1, 0, 0], [0, 1, 0], hr3, finv2);

disp(' ');
disp('***** Estimated VARMA(1,0)(0,1) Model using the HR method  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'phi';
strt = 1;
mprintar(strv.phis3(:, :, 2), in, tit, strt);
disp(' ');
tit = 'th';
strt = 12;
mprintar(strv.thetas3(:, :, 13), in, tit, strt);
disp('press any key to continue')
pause


%estimate using exact ML
%setup model
Phi = eye(2);
Th(:, :, 1) = eye(2);
Th(:, :, 2) = strv.thetas3(:, :, 13);
phi = strv.phis3(:, :, 1:2);
th = eye(2);
Sigma = strv.sigmar3;

[str, ferror] = suvarmapqPQ(phi, th, Phi, Th, Sigma, freq);

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
tit = 'phi';
strt = 1;
mprintar(phif(:, :, 2), in, tit, strt);
disp(' ')
tit = 'th';
strt = 12;
mprintar(Thf(:, :, 2), in, tit, strt);
disp(' ')
tit = 'Sigma';
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
tit = 'tv-phi';
strt = 1;
mprintar(phitvf(:, :, 2), in, tit, strt);
disp(' ')
tit = 'tv-th';
strt = 12;
mprintar(Thtvf(:, :, 2), in, tit, strt);
disp(' ')
tit = 'tv-L';
mprintar(Ltvf, in, tit);
