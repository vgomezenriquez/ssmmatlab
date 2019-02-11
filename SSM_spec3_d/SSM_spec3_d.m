%*****************************************************************
%
%    SPECTRAL ANALYSIS: US INDUSTRIAL PRODUCTION, CONSUMPTION AND HOURS
%
% Series: Cycles of the US IPI, consumption and working hours (monthly data)
% Time span: 1953.M4 - 2007.M9
%*****************************************************************

clear;
clc;

y = load(fullfile('data', 'USc3.dat'));
% Settings for spectral analysis:
per = 12; %number of seasons
win = 2; %Parzen window
corlag = 40; %number of correlation and cross correlations
graph = 1; %plot results
vnames = {'US IPI', 'US CONSUMPTION', 'US HOURS'}; %names for the series

spr = spectralan(y, per, win, corlag, graph, vnames);
