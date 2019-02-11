[<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/banner.png" width="888" alt="Visit QuantNet">](http://quantlet.de/)

## [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/qloqo.png" alt="Visit QuantNet">](http://quantlet.de/) **SSM_CovFac_d** [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/QN2.png" width="60" alt="Visit QuantNet 2.0">](http://quantlet.de/)

```yaml

Name of QuantLet: SSM_CovFac_d

Published in: Linear Time Series With MATLAB and Octave

Description: 'The autocovariances of a vector moving averrage process of order one are first computed and
               then the spectral factorization of these covariances is performed using two methods.'

Keywords: time-series, spectral factorization, vector moving average, autocovariances, polynomial methods, state space methods

Author: Víctor Gómez

Submitted: Fri, January 25 2019 by Víctor Gómez

```

### MATLAB Code
```matlab

%
%script file to illustrate two covariance factorization methods. The first
%one is Tunnicliffe-Wilson method and the second one uses the DARE.
%
clear

phi(:, :, 1) = eye(2);
th(:, :, 1) = eye(2);
th(:, :, 2) = [6.5, -2.; 15, -4.5];
Sigma = [4., 1.; 1., 1.];

nc = 2;
[c, ierror] = macgf(phi, th, Sigma, nc);
disp('Autocovariance matrices of lags 0 to 1:')
for i = 1:2
    disp(c(:, :, i))
end

disp('covariance factorization using Tunnicliffe-Wilson method')

[Omega, Theta, ierror, iter, normdif] = pmspectfac(c, 50)

disp('press any key to continue')
pause

disp('compute autocovariances of the invertible model')
phin(:, :, 1) = eye(2);
thn(:, :, 1) = eye(2);
thn(:, :, 2) = Theta;
Sigman = Omega;

nc = 2;
[cn, ierror] = macgf(phin, thn, Sigman, nc);
disp('Autocovariance matrices of lags 0 to 1:')
for i = 1:2
    disp(cn(:, :, i))
end
disp('press any key to continue')
pause

disp('covariance factorization using the DARE')
[Thetad, Omegad] = ssmspectfac(c)

```

automatically created on 2019-02-11