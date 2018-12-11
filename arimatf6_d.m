%Example of automatic model identification, estimation and forecasting of a
%transfer function model.
%The model is a simulated trasfer function model (SCA)
% TX2 AND TY2 ARE TWO SIMULATED SERIES.  THE MODEL IS
%      (1 - B)Y(T) = (3.0B - 2.0B**2)(1-B)X(T) + (1 - 0.7B)A(T)
% THERE ARE 130 OBSERVATIONS IN THIS DATA SET.
%
%The specification file, tf2.m, is in the subdirectory spec. In this
%specification file, the instructions for the program are given. The
%default values for the program are in the script file arimadefval.m.
%

out = arimaestos('tf2');

disp(' ')
disp('Details on the identified model, its estimation, and the forecasts')
disp('are in the file "tf2.txt" in the subdirectory "results"')
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