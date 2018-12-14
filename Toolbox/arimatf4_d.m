%Example of automatic model identification and estimation of a transfer
%function model.
%The model is that of Box and Jenkins, (1976)
% SALES DATA WITH LEADING INDICATOR
% THERE ARE 150 OBSERVATIONS
% The identified and estimated model is (SCA Liu, 2005, p. 5.14):
% (1-B)y_t = 0.035 + 4.726*B^3/(1-0.724*B)(1-B)x_t + (1-0.626*B)a_t.
%
%The specification file, bjsales.m, is in the subdirectory spec. In this
%specification file, the instructions for the program are given. The
%default values for the program are in the script file arimadefval.m.
%

out = arimaestos('bjsales');

disp(' ')
disp('Details on the identified model, its estimation, and the forecasts')
disp('are in the file "bjsales.txt" in the subdirectory "results"')
disp('press any key to continue')
pause

disp(' ')
disp('Structure out contains the following fields:')
disp('.title: output series name')
disp('.nziyip: number of observations, initial year and initial month')
disp('  or quarter');
disp('.freq: series frequency')
disp('.orig: original series')
out
disp('press any key to continue')
pause

disp('.model: contains fields with information about the identified and ')
disp('  estimated ARIMA model when finite linear approximations to the ')
disp('  input filters are used:')
out.model
disp('press any key to continue')
pause

disp('.tfmodel: contains information about the identified and ')
disp(' estimated transfer function model and the forecasts:')
out.tfmodel