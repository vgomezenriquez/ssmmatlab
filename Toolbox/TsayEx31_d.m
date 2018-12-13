%script file for example3.1 in Tsay (2014)
%

data = load(fullfile('data', 'm-dec15678-6111.dat'));
x = log(data(:, 2:6)+1) * 100;
% size(x)
rtn = x(:, [2, 5]);
tdx = [1:612] / 12 + 1961;

plot(tdx, rtn(:, 1))
xlabel('year');
ylabel('d5');
pause

plot(tdx, rtn(:, 2))
xlabel('year');
ylabel('d8');
pause

close all

lag = 6;
ic = 1;
str = mautcov(rtn, lag, ic);
disp('Correlation matrix at lag 0:')
disp(str.r0)
