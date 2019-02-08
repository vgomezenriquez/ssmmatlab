% script file to first simulate a normally distributed univariate white
% noise series and then compute some sample covariances and Ljung-Box
% statistics

clear

l = 50;    %number of initial observations to be discarded in the
           %simulated series
N = 150;   %number of observations of the simulated series
s = 1;     %number of outputs
seed = 20; %this is to always generate the same series

%polynomial matrices Phi and Theta
Phi(:, :, 1) = eye(s);
Theta(:, :, 1) = eye(s);
%we add the following term to pretend we simulate a VARMA model, but
%in fact we are simulating white noise
Phi(:, :, 2) = zeros(s);

%covariance matrix of the a_t innovartions
S = 2.;

%simulate v_t
[v, ferror] = varmasim(l, N, Phi, Theta, S, seed);

%sample autocovariances and autocorrelations
lag = 10;
ic = 1;
nr = 0;
[c0, cv, r] = autcov(v, lag, ic);
%Q statistics
nv = size(v, 1);
orders = 1:lag;
[qstat, pval] = lbs(nv, orders, r, nr);
disp('Q statistics:')
disp(qstat)

disp('p-values of Q statistics:')
disp(pval)
