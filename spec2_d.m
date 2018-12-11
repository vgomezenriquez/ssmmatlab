%*****************************************************************
%
%    SPECTRAL ANALYSIS: GERMAN INDUSTRIAL PRODUCTION
%
% Series: German IPI cycle (quarterly data)
% Time span: 1970.Q1 - 2011.Q3
%*****************************************************************

clear;
clc;

y = load(fullfile('data', 'PRGer.dat'));
% Settings for spectral analysis:
per = 4; %number of seasons
win = 2; % Parzen window
corlag = 30;
graph = 1;
vnames = {'German IPI'};

spr = spectralan(y, per, win, corlag, graph, vnames);
