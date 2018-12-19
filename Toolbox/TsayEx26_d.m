%script file for example 2.6 in Tsay (2014)
%

data = load(fullfile('data', 'q-gdp-ukcaus.dat'));
gdp = log(data(:, 3:5));
% size(gdp)

z = diferm(gdp, 1); % Growth rate
z = z * 100; % Percentage growth rate

disp(' ')
disp('estimation of a VAR(2) model for the series')
disp('')
disp('press any key to continue')
pause

%estimate VAR(2)
nlag = 2;
res = var_est(z, nlag);

nvar = 3;
lag = 12;
ic = 1;
nr = nlag * nvar^2;
disp('******** VAR residuals:     ********');
str = mautcov(res.resid, lag, ic, nr);
disp('Q statistics:')
disp(str.qstat)
disp('press any key to continue')
pause
disp('p-values of Q statistics:')
disp(str.pval)
[m, n] = size(str.pval);
t = 1:m;
plot(t, str.pval, t, 0.05*ones(1, m))
legend('p-values of Q statistics:')
disp('press any key to continue')
pause
close all
