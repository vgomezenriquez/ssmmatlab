[<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/banner.png" width="888" alt="Visit QuantNet">](http://quantlet.de/)

## [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/qloqo.png" alt="Visit QuantNet">](http://quantlet.de/) **SSM_TsayEx37_d** [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/QN2.png" width="60" alt="Visit QuantNet 2.0">](http://quantlet.de/)

```yaml

Name of QuantLet: SSM_TsayEx37_d

Published in: Linear Time Series With MATLAB and Octave

Description: 'The first three covariance matrices of a VARMA(2,1) are computed. Example 3.7 in Tsay (2014).'

Keywords: time-series, VARMA model, covariance matrices, polynomial matrices, cross-correlation matrices

Author: Víctor Gómez

Submitted: Tue, December 18 2018 by Víctor Gómez

```

### MATLAB Code
```matlab

%script file for example3.7 in Tsay (2014)
%

phi(:, :, 1) = eye(2);
phi(:, :, 2) = -[0.816, -0.623; -1.116, 1.074];
phi(:, :, 3) = -[-0.643, 0.592; 0.615, -0.133];
th(:, :, 1) = eye(2);
th(:, :, 2) = -[0, -1.248; -0.801, 0];
Sigma = [4, 2; 2, 5]


nc = 3;
[c, ierror] = macgf(phi, th, Sigma, nc);
disp('Autocovariance matrices of lags 0 to 2:')
for i = 1:3
    disp(c(:, :, i))
end

```

automatically created on 2019-02-11