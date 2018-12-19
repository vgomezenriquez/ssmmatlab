%script file for example3.4 in Tsay (2014)
%

data = load(fullfile('data', 'm-ibmko-0111.dat'));
x = log(data(:, 2:3)+1.) * 100;
[nx, mx] = size(x);
rtn = x;
tdx = [1:nx] / 12 + 2001;

subplot(2, 1, 1)
plot(tdx, rtn(:, 1))
xlabel('year');
ylabel('ln-rtn');
axis('tight');
subplot(2, 1, 2)
plot(tdx, rtn(:, 2))
xlabel('year');
ylabel('ln-rtn');
axis('tight');
disp('press any key to continue')
pause
close all


%ccm matrices
lag = 10;
ic = 1;
nr = 0;
str = mautcov(rtn, lag, ic, nr);
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

%difference data
y = diferm(rtn, 1);

%estimate an MA(1) model using the Hannan-Rissanen method
disp(' ')
disp('estimate an MA(1) model using the Hannan-Rissanen method: ')
disp('press any key to continue')
pause

freq = 1;
xx = [];
hr3 = 0;
finv2 = 1;
[strv, ferror] = estvarmaxpqrPQR(y, xx, freq, [0, 1, 0], [0, 0, 0], hr3, finv2);

disp(' ');
disp('***** Estimated  MA(1)  Model  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'mu:';
mprintar(strv.mus3', in, tit);
disp('press any key to continue')
pause
tit = 'th:';
strt = 1;
mprintar(strv.thetas3(:, :, 2), in, tit, strt);
disp('press any key to continue')
pause

tit = 'tv-th:';
strt = 1;
mprintar(strv.thetatv3(:, :, 2), in, tit, strt);
disp('press any key to continue')
pause
tit = 'Sigma:';
mprintar(strv.sigmar3, in, tit);
disp('press any key to continue')
pause

disp(' ');
disp('***** Estimate Model using the conditional method  *****');
disp(' ');
disp('press any key to continue')
pause

%estimate the model using the conditional method
[xvfc, strc, ferrorc] = mconestim(y, xx, strv);

disp(' ');
disp('***** Estimated Model using the conditional method  *****');
disp(' ');
clear in
in.fid = 1;
tit = 'th:';
strt = 1;
mprintar(strc.thetascon(:, :, 2), in, tit, strt);

disp(' ');
tit = 'Sigma:';
mprintar(strc.sigmarcon, in, tit);

disp(' ')
disp('t-values: ')
in.fmt = char('%12.4f');
disp(' ');
tit = 'tv-th:';
strt = 1;
mprintar(strc.thetatvcon(:, :, 2), in, tit, strt);

disp('press any key to continue')
pause


%estimate using the exact method
disp(' ')
disp('estimation using the exact method')
disp('press any key to continue')
pause
Y = [];
[xvf, strx, ferrorc] = mexactestim(y, xx, strv, Y);

disp(' ');
disp('***** Estimated Model using the exact method  *****');
disp(' ');
clear in
in.fid = 1;
tit = 'th:';
strt = 1;
mprintar(strx.thetasexct(:, :, 2), in, tit, strt);

disp(' ');
tit = 'Sigma:';
mprintar(strx.sigmarexct, in, tit);

disp(' ')
disp('t-values: ')
in.fmt = char('%12.4f');
disp(' ');
tit = 'tv-th:';
strt = 1;
mprintar(strx.thetatvexct(:, :, 2), in, tit, strt);

disp('press any key to continue')
pause


disp('eigenvalues of both the conditional and the exact method')
eig(-strc.thetascon(:, :, 2))
eig(-strx.thetasexct(:, :, 2))
