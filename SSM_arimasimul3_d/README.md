[<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/banner.png" width="888" alt="Visit QuantNet">](http://quantlet.de/)

## [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/qloqo.png" alt="Visit QuantNet">](http://quantlet.de/) **SSM_arimasimul3_d** [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/QN2.png" width="60" alt="Visit QuantNet 2.0">](http://quantlet.de/)

```yaml

Name of QuantLet: SSM_arimasimul3_d

Published in: Linear Time Series With MATLAB and Octave

Description: 'Firstly, a normally distributed white noise is simulated. Then, the Ljung-Box statistic is computed for several lags.'

Keywords: time-series, white noise, simulation, Ljung-Box statistic, simulation, p-values

Author: Víctor Gómez

Submitted: Wed, December 19 2018 by Víctor Gómez

```

### MATLAB Code
```matlab

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

```

automatically created on 2019-02-11