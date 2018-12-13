%script file for example 2.7 in Tsay (2014)
%

data = load(fullfile('data', 'q-gdp-ukcaus.dat'));
gdp = log(data(:, 3:5));
% size(gdp)

z = diferm(gdp, 1); % Growth rate
z = z * 100; % Percentage growth ra

disp(' ')
disp('estimation of a VAR(2) model for the series')
disp('')
disp('press any key to continue')
pause

%estimate VAR(2)
nlag = 2;
res = var_est(z, nlag);
disp(' ');
disp('***** Estimated VAR Model  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'AR';
strt = 1;
mprintar(res.phi(:, :, 2:3), in, tit, strt);
disp(' ')
tit = 'Constant';
mprintar(res.const', in, tit);

disp(' ');
disp('*****t-values  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'tv-AR';
strt = 1;
mprintar(res.phitv(:, :, 2:3), in, tit, strt);
disp(' ')
tit = 'tv-Constant';
mprintar(res.consttv', in, tit);
disp(' ');
tit = 'Sigma:';
mprintar(res.sigmar, in, tit);
disp('press any key to continue')
pause

%eliminate some nonsignificant  parameters in the model using function
%estvarmaxpqrPQR.
disp(' ')
disp('Elimination of nonsignificant parameters in the VAR(2) model ')
disp('')
disp('press any key to continue')
pause

freq = 1;
xx = [];
hr3 = 1;
finv2 = 0;
nsig = [1, 0];
tsig = [1.96, 0.];
mstainv = 0;
[strvr, ferror] = estvarmaxpqrPQR(z, xx, freq, [2, 0, 0], [0, 0, 0], hr3, finv2, ...
    mstainv, nsig, tsig);
disp(' ');
disp('***** Estimated  simplified VAR Model  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'AR';
strt = 1;
mprintar(strvr.phis(:, :, 2:3), in, tit, strt);
disp(' ')
tit = 'Constant';
mprintar(strvr.mu', in, tit);

disp(' ');
disp('*****t-values  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'tv-AR';
strt = 1;
mprintar(strvr.phitv(:, :, 2:3), in, tit, strt);
disp(' ')
tit = 'tv-Constant';
mprintar(strvr.mutv', in, tit);
disp(' ');
tit = 'Sigma:';
mprintar(strvr.sigmar2, in, tit);
disp('press any key to continue')
pause


nvar = 3;
lag = 12;
ic = 1;
nr = strvr.nparm;
disp('******** VAR residuals:     ********');
str = mautcov(strvr.resid2, lag, ic, nr);
disp('Q statistics:')
disp(str.qstat)
disp('press any key to continue')
pause
disp('p-values of Q statistics:')
disp(str.pval)
[m, n] = size(str.pval);
t = 1:m;
plot(t, 0.05*ones(1, m), t, str.pval)
legend('p-values of Q statistics:')
disp('press any key to continue')
pause
close all
