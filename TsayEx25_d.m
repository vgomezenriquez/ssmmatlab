%script file for example 2.5 in Tsay (2014)
%

data = load(fullfile('data', 'q-gdp-ukcaus.dat'));
gdp = log(data(:, 3:5));
% size(gdp)

y = diferm(gdp, 1); % Growth rate

%VAR order identification
disp(' ');
disp('***** VAR order identification  *****');
disp(' ');
prt = 1;
minlag = 0;
maxlag = 13;
lagsopt = varident(y, maxlag, minlag, prt);
