%Example of estimation of the flour price series using a
%VARMA(1,1)  model.
%
% Monthly Flour Price Indices for Three U.S. cities, Buffalo, Minneapolis,
% and Kansas City, for the Period August 1972 Through November 1980. The
% three series are I(1) and it seems that there are no cointegration
% relationships. These series have been used by Tiao y Tsay (1989) and
% Lütkepohl and Poskitt (1996).

clear
%load data
y = load(fullfile('data', 'flour-price.dat'));
x = [];
%logs are taken
y = log(y);

lag = 20;
cw = 1.96;
ds = 0;
dr = 1;
seas = 1;
tname = {'flour price 1', 'flour price 2', 'flour price 3'};
for i = 1:3
    c0 = sacspacdif(y(:, i), tname(i), dr, ds, seas, lag, cw);
    pause
end
close all


freq = 1;
xx = [];
hr3 = 0;
finv2 = 1;
[strv, ferror] = estvarmaxpqrPQR(y, xx, freq, [1, 1, 0], [0, 0, 0], hr3, finv2);

%estimate using the conditional method
disp(' ')
disp('estimation using the conditional method')
disp('model is VARMA(1,1)')
disp('press any key to continue')
pause
[xvfc, strc, ferrorc] = mconestim(y, xx, strv);

disp(' ');
disp('***** Estimated Model using the conditional method  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'mu:';
mprintar(strc.muscon', in, tit);

tit = 'phi';
strt = 1;
mprintar(strc.phiscon(:, :, 2), in, tit, strt);

tit = 'th';
strt = 1;
mprintar(strc.thetascon(:, :, 2), in, tit, strt);

disp(' ');
tit = 'Sigma:';
mprintar(strc.sigmarcon, in, tit);

disp(' ')
disp('t-values: ')
disp(' ');
tit = 'tv-mu:';
mprintar(strc.mutvcon', in, tit);

tit = 'tv-phi';
strt = 1;
mprintar(strc.phitvcon(:, :, 2), in, tit, strt);

tit = 'tv-th';
strt = 1;
mprintar(strc.thetatvcon(:, :, 2), in, tit, strt);
disp('press any key to continue')
pause


lag = 24;
ic = 1;
nr = strc.nparm;
disp(' ')
disp('******** Conditional Residuals:     ********');
str = mautcov(strc.residcon, lag, ic, nr);
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

%estimate using the exact method
disp(' ')
disp('estimation using the exact method')
disp('model is VARMA(1,1)')
disp('press any key to continue')
pause
[ny, s] = size(y);
Y = eye(s);
[xvf, strx, ferrorc] = mexactestim(y, xx, strv, Y);


disp(' ');
disp('***** Estimated Model using the exact method  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'mu:';
mprintar(strx.musexct', in, tit);

tit = 'phi';
strt = 1;
mprintar(strx.phisexct(:, :, 2), in, tit, strt);

tit = 'th';
strt = 1;
mprintar(strx.thetasexct(:, :, 2), in, tit, strt);

disp(' ');
tit = 'Sigma:';
mprintar(strx.sigmarexct, in, tit);

disp(' ')
disp('t-values: ')
disp(' ');
tit = 'tv-mu:';
mprintar(strx.mutvexct', in, tit);

tit = 'tv-phi';
strt = 1;
mprintar(strx.phitvexct(:, :, 2), in, tit, strt);

tit = 'tv-th';
strt = 1;
mprintar(strx.thetatvexct(:, :, 2), in, tit, strt);
disp('press any key to continue')
pause
