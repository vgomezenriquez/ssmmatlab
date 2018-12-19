%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example  6.2 (continued)  in Tsay (2014), pp. 353-357
%
%Series are:  1) monthly U.S. regular conventional gas price z1t and
% 2)  heating oil price z2t of New York Harbor. Both series are measured in dollars
% per gallon. These prices depend on the crude oil and natural gas prices;  3) x1t
% the spot oil price of West Texas Intermediate, dollars per barrel, and  4) x2t the
% natural gas price of Henry Hub, LA, measured in dollars per million BTU. Thus,
% yt =(z1t, z2t), xt = (x1t, x2t), and k = m = 2. The sample period is from
% November 1993 to August 2012.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

da = load(fullfile('data', 'm-gasoil.dat'));
yt = da(:, 3:6);
tdx = da(:, 1) + da(:, 2) / 12;


subplot(2, 2, 1)
plot(tdx, yt(:, 1))
xlabel('time');
ylabel('Greg');
axis('tight');
subplot(2, 2, 2)
plot(tdx, yt(:, 3))
xlabel('time');
ylabel('oilp');
axis('tight');
subplot(2, 2, 3)
plot(tdx, yt(:, 2))
xlabel('time');
ylabel('hoil');
axis('tight');
subplot(2, 2, 4)
plot(tdx, yt(:, 4))
xlabel('time');
ylabel('gasp');
axis('tight');
disp('press any key to continue')
pause
close all

zt = yt(:, 1:2);
xt = yt(:, 3:4);

%perform multivariate linear regression
nlag = 0;
test = 0;
strr = var_est(zt, nlag, test, xt);

disp(' ');
disp('***** Estimated Regression Model  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'Beta';
strt = 1;
mprintar(strr.betavar(1:2, :)', in, tit, strt);
disp(' ')
tit = 'Constant';
mprintar(strr.const', in, tit);
disp(' ')
tit = 'tv-Beta';
strt = 1;
mprintar(strr.tvvar(1:2, :)', in, tit, strt);
disp(' ')
tit = 'tv-Constant';
mprintar(strr.consttv', in, tit);
disp('press any key to continue')
pause

% residuals without mean
strr.betavar(2, 1) = 0.; %eliminate an insignificant parameter
res = zt - xt * strr.betavar(1:2, :);

%identify a VAR for the residuals
disp('identify a VAR model for the residuals')
disp(' ')
prt = 1;
minlag = 0;
maxlag = 13;
lagsopt = varident(res, maxlag, minlag, prt);
disp('press any key to continue')
pause


disp('estimate a simplified VAR(2) model for the residuals using the HR method')
disp(' ')
disp('press any key to continue')
pause
%estimate a simplified VAR(2) model for the residuals
freq = 1;
hr3 = 0;
finv2 = 1;
mstainv = 1;
nsig = [1, 0];
tsig = [1.2, 1.];
[strvs, ferror] = estvarmaxpqrPQR(res, xt, freq, [2, 0, 0], [0, 0, 0], hr3, finv2, mstainv, nsig, tsig);

disp(' ');
disp('***** Estimated VAR Model  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'AR';
strt = 1;
mprintar(strvs.phis3(:, :, 2:3), in, tit, strt);
disp(' ')
tit = 'Constant';
mprintar(strvs.mus3', in, tit);
disp(' ')
tit = 'tv-AR';
strt = 1;
mprintar(strvs.phitv3(:, :, 2:3), in, tit, strt);
disp(' ')
tit = 'tv-Constant';
mprintar(strvs.mutv3', in, tit);
disp(' ')
tit = 'Sigma';
mprintar(strvs.sigmar3, in, tit);
disp('press any key to continue')
pause


%compute autocovariance and autocorrelation matrices of  residuals
lag = 24;
ic = 1;
nr = 0; %nr=strx.nparm;
disp(' ')
disp('******** Residuals:     ********');
stre = mautcov(strvs.resid3, lag, ic, nr);
disp('Correlation matrix at lag 0:')
disp(stre.r0)
disp('Q statistics:')
disp(stre.qstat)

disp('p-values of Q statistics:')
disp(stre.pval)
[m, n] = size(stre.pval);
t = 1:m;
plot(t, stre.pval, t, 0.05*ones(1, m))
legend('p-values of Q statistics:')
pause
close all

disp(' ')
disp('estimation using the exact method')
disp('press any key to continue')
pause
%estimate using exact ML
%setup model
pr = 2;
qr = 0;
Phi = eye(2);
Th(:, :, 1) = eye(2);
th(:, :, 1) = eye(2);
phi = strvs.phis3(:, :, 1:pr+1);
Sigma = strvs.sigmar3;

[str, ferror] = suvarmapqPQ(phi, th, Phi, Th, Sigma, freq);

%eliminate insignificant parameters
mx = 2;
for k = 2:pr + 1
    for i = 1:mx
        for j = 1:mx
            if phi(i, j, k) == 0
                str.phin(i, j, k) = 0;
            end
        end
    end
end
[str, ferror] = fixvarmapqPQ(str);

%regression matrix
Y = [repmat(eye(2), size(zt, 1), 1), kron(xt, eye(2))];


%estimate model using the exact method
result = varmapqPQestim(zt, str, Y);

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
mprintar(phif(:, :, 2:pr+1), in, tit, strt);
disp(' ')
tit = 'Beta';
Beta = [result.h(3:4, :), result.h(5:6, :)];
mprintar(Beta, in, tit);
disp(' ')
tit = 'Constant';
ct = result.h(1:2)';
mprintar(ct, in, tit);
disp(' ')
tit = 'tv-phi';
strt = 1;
mprintar(phitvf(:, :, 2:pr+1), in, tit, strt);
disp(' ')
tit = 'tv-Beta';
Beta = [result.tvr(3:4, :), result.tvr(5:6, :)];
mprintar(Beta, in, tit);
disp(' ')
tit = 'tv-Constant';
ct = result.tvr(1:2)';
mprintar(ct, in, tit);
disp(' ')
tit = 'Sigma:';
mprintar(Sigmar, in, tit);
