%script file for example3.3 in Tsay (2014)
%

data = load(fullfile('data', 'm-dec15678-6111.dat'));
x = log(data(:, 2:6)+1) * 100;
% size(x)
rtn = x(:, [2, 5]);
tdx = [1:612] / 12 + 1961;

subplot(2, 1, 1)
plot(tdx, rtn(:, 1))
xlabel('year');
ylabel('d5');
axis('tight');
subplot(2, 1, 2)
plot(tdx, rtn(:, 2))
xlabel('year');
ylabel('d8');
axis('tight');
disp('press any key to continue')
pause
close all


%estimate an MA(1) model using the Hannan-Rissanen method
disp(' ')
disp('estimate an MA(1) model using the Hannan-Rissanen method: ')
disp('press any key to continue')
pause

freq = 1;
xx = [];
hr3 = 0;
finv2 = 1;
[strv, ferror] = estvarmaxpqrPQR(rtn, xx, freq, [0, 1, 0], [0, 0, 0], hr3, finv2);

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
disp('***** Estimate Model using the exact method  *****');
disp(' ');
disp('press any key to continue')
pause

%estimate the model using the exact method
[ny, s] = size(rtn);
Y = repmat(eye(s), ny, 1);
[xvf, strx, ferror] = mexactestim(rtn, xx, strv, Y);

disp(' ');
disp('***** Estimated Model using the exact method  *****');
disp(' ');
clear in
in.fid = 1;
tit = 'mu:';
mprintar(strx.musexct', in, tit);
disp('press any key to continue')
pause
tit = 'th:';
strt = 1;
mprintar(strx.thetasexct(:, :, 2), in, tit, strt);

disp(' ');
tit = 'Sigma:';
mprintar(strx.sigmarexct, in, tit);

disp(' ')
disp('t-values: ')
in.fmt = char('%12.4f');
tit = 'tv-mu:';
mprintar(strx.mutvexct', in, tit);
disp('press any key to continue')
pause

disp(' ');
tit = 'tv-th:';
strt = 1;
mprintar(strx.thetatvexct(:, :, 2), in, tit, strt);

disp('press any key to continue')
pause


% disp(' ')
% disp('******** Computation of Recursive Residuals   ********');
%compute recursive residuals
%set up regression matrices
X = eye(s);
W = [];
Sigmax = strx.sigmarexct;
[L, p] = chol(Sigmax, 'lower');
Lm1 = pinv(L);
%set up system matrices
T = strx.Fsexct;
Z = strx.Hsexct;
G = Lm1;
H = strx.Ksexct * Lm1;
%set up initial conditions
ndelta = 0; %number of unit roots
[ins, ii, ferror] = incossm(T, H, ndelta);
disp(' ')
disp('******** Computation of OLS Residuals   ********');
%compute OLS residuals
g = strx.musexct;
[KKP, PT, recrs, recr, srecr, t1, A1, LP1, KG] = scakfffsqrt(rtn, X, Z, G, W, T, H, ins, ii, g);
%plot OLS residuals
plot(recr(:, 1)), legend('olsres(:,1)'), pause
plot(recr(:, 2)), legend('olsres(:,2)'), pause
close all
%compute autocovariance and autocorrelation matrices of rec. residuals
lag = 24;
ic = 1;
nr = strx.nparm;
disp(' ')
disp('******** OLS Residuals:     ********');
str = mautcov(recr, lag, ic, nr);
disp('Correlation matrix at lag 0:')
disp(str.r0)
disp('Q statistics:')
disp(str.qstat)

disp('p-values of Q statistics:')
disp(str.pval)
[m, n] = size(str.pval);
t = 1:m;
plot(t, 0.05*ones(1, m), t, str.pval)
legend('p-values of Q statistics:')
disp('press any key to continue')
pause
close all
