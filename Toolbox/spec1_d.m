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
