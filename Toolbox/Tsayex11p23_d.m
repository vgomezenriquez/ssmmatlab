%script file to carry out exercise 1.1 on page 23 of Tsay 2014 book
%regarding the simulation of a multivariate series following a VAR(1) model.
%

clear
l = 50; %number of initial observations to be discarded in the simulated series
N = 300; %number of observations of the simulated series
s = 2; %number of variables
% rng default   %default startup generator:
seed = 20; %this is to always generate the same series

phi(:, :, 1) = eye(s);
th(:, :, 1) = eye(s);

%VAR(1) matrix
phi(:, :, 2) = -[0.8, 0.4; -0.3, 0.6];
%covariance matrix of the innovartions
S = [2.0, 0.5; 0.5, 1.0];

%simulate VAR(1)
[z, ferror] = varmasim(l, N, phi, th, S, seed);

nvar = 2;
lag = 10;
ic = 1;
nr = 0;
disp('******** CCM matrices:     ********');
str = mautcov(z, lag, ic, nr);
disp('Q statistics:')
disp(str.qstat)
disp('p-values of Q statistics:')
disp(str.pval)
disp('press any key to continue')
pause

disp(' ');
disp('***** VAR order identification  *****');
disp(' ');
maxlag = 6;
minlag = 1;
prt = 1;
lagsopt = varident(z, maxlag, minlag, prt);
nlag = 1;
res = var_est(z, nlag);
disp('press any key to continue')
pause


disp(' ');
disp('***** Estimated VAR Model  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'AR';
strt = 1;
mprintar(res.phi(:, :, 2), in, tit, strt);
disp(' ')
tit = 'Constant';
mprintar(res.const', in, tit);


disp(' ');
disp('***** Estimated t-values  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'tv-AR';
strt = 1;
mprintar(res.phitv(:, :, 2), in, tit, strt);
disp(' ')
tit = 'tv-Constant';
mprintar(res.consttv', in, tit);
disp(' ');
tit = 'Sigma:';
mprintar(res.sigmar, in, tit);
