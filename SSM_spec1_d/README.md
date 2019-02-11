[<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/banner.png" width="888" alt="Visit QuantNet">](http://quantlet.de/)

## [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/qloqo.png" alt="Visit QuantNet">](http://quantlet.de/) **SSM_spec1_d** [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/QN2.png" width="60" alt="Visit QuantNet 2.0">](http://quantlet.de/)

```yaml

Name of QuantLet: SSM_spec1_d

Published in: Linear Time Series With MATLAB and Octave

Description: 'The smoothed periodogram of the cycle of the quarterly German IPI series is computed
               and displayed.'

Keywords: time-series, smoothed periodogram, cycle, ARIMA model, canonical decomposition

Author: Víctor Gómez

Submitted: Fri, January 25 2019 by Víctor Gómez

```

![Picture1](spgipi.png)

### MATLAB Code
```matlab

%*****************************************************************
%
%    SPECTRAL ANALYSIS: GERMAN INDUSTRIAL PRODUCTION
%
% Series: German IPI cycle (quarterly data)
% Time span: 1970.Q1 - 2011.Q3
%*****************************************************************

clear;
clc;

x = load(fullfile('data', 'PRGer.dat'));
rfname = 'PRGer';

% Settings for spectral analysis:
win = 2; % Parzen window

% compute smothed periodogram of IPI
[spx, frq] = periodg(x, win);

%plot spectrum
per = 4; % data frequency
cc = ones(100, 2); %lines corresponding to cycle frequency band
% frequency band for which the results are displayed;
% it corresponds to business cycle periodicities (periods between 1.5 and 8
% years)
cc(:, 1) = cc(:, 1) * ((2 * pi) / (per * 1.5));
cc(:, 2) = cc(:, 2) * ((2 * pi) / (per * 8));
ll = (max(spx) - min(spx)) / 99;
dd = min(spx):ll:max(spx);

figure
plot(frq, spx, cc(:, 1), dd, cc(:, 2), dd)
legend(['smoothed spectrum ', rfname])
disp('strike a key to continue')
pause
% Close figure
close all

```

automatically created on 2019-02-11