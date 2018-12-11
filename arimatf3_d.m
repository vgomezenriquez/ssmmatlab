%Example of automatic model identification and estimation of an ARIMA
%model. The program performs automatic outlier detection.
% Series is Variety stores sales (N=153), used by Hillmer, Bell and Tiao
% (1983)
%
%The specification file, vstores.m, is in the subdirectory spec. In this
%specification file, the instructions for the program are given. The
%default values for the program are in the script file arimadefval.m.
%

out = arimaestos('vstores');

disp(' ')
disp('Details on the identified model, its estimation, and outlier ')
disp('detection are in the file "vstores.txt" in the subdirectory')
disp('"results"')
disp('press any key to continue')
pause

disp('Structure out contains the following fields:')
disp('.title: series name')
disp('.nziyip: number of observations, initial year and initial month')
disp('  or quarter');
disp('.freq: series frequency')
disp('.orig: original series')
out
disp('press any key to continue')
pause

disp('.model: contains fields with information about the identified and ')
disp('  estimated model and the forecasts:')
out.model
disp('.tfmodel: contains information about transfer function models ')