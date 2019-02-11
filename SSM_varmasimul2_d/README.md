[<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/banner.png" width="888" alt="Visit QuantNet">](http://quantlet.de/)

## [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/qloqo.png" alt="Visit QuantNet">](http://quantlet.de/) **SSM_varmasimul2_d** [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/QN2.png" width="60" alt="Visit QuantNet 2.0">](http://quantlet.de/)

```yaml

Name of QuantLet: SSM_varmasimul2_d

Published in: Linear Time Series With MATLAB and Octave

Description: 'Firstly, a time series following a MA(1) model is simulated. Then, some cross correlation matrices are computed.'

Keywords: time-series, VARMA model, simulation, cross correlation matrices, polynomial matrices

Author: Víctor Gómez

Submitted: Wed, December 19 2018 by Víctor Gómez


```

### MATLAB Code
```matlab

%script file to simulate a VARMA model
%
% Let the VARMA model be
%
%   Phi(B) y_t = Theta(B) a_t,
%
% where B is the backshift operator, By_t = y_{t-1}.
%
clear 

l = 50;     %number of initial observations to be discarded in the
            %simulated series
N = 300;    %number of observations of the simulated series
s = 3;      %number of outputs
seed = 20;  %this is to always generate the same series

%polynomial matrices Phi and Theta
Phi(:, :, 1) = eye(s);
Theta(:, :, 1) = eye(s);
Theta(:, :, 2) = -[0.8, 0.4, 0.; -0.3, 0.6, 0.; 0., 0., -1.];

%covariance matrix of the a_t innovartions
S = [2.0, 0.5, 0.; 0.5, 1.0, 0.; 0., 0., .3];

%simulate v_t
[v, ferror] = varmasim(l, N, Phi, Theta, S, seed);

disp('cross correlation matrices')
%ccm matrices
lag = 6;
ic = 1;
str = mautcov(v, lag, ic);
disp('Correlation matrix at lag 0:')
disp(str.r0)

```

automatically created on 2019-02-11