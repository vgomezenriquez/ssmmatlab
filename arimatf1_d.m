%Example of automatic model identification, estimation and forecasting of
%an ARIMA model.
%Series is airline series from Box and Jenkins (1976)
%
%The specification file, bjsgairl.m, is in the subdirectory spec. In this
%specification file, the instructions for the program are given. The
%default values for the program are in the script file arimadefval.m.
%

out = arimaestos('bjsgairl');

disp(' ')
disp('Details on the identified model, its estimation, and the forecasts')
disp('are in the file "bjsgairl.txt" in the subdirectory "results"')
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
