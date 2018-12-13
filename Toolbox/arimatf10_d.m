%Example of an ARIMA model with complex seasonal patterns.
%Series is Turkey Electricity Demand, data analyzed in De Livera,
%Hyndman and Snyder (2013), Journal of the American Statistical
%Association, ''Forecasting Time Series With Complex Seasonal Patterns
%Using Exponential Smoothing'', 106, 1513-1527.
%
%The specification file, turkey_elec.m, is in the subdirectory spec. In
%this specification file, the instructions for the program are given. The
%default values for the program are in the script file arimadefval.m.
%

out = arimaestos('turkey_elec');

disp(' ')
disp('Details on the identified and estimated model are in the')
disp('file "turkey_elec.txt" in the subdirectory "results"')
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
disp('  estimated model:')
out.model
disp('.tfmodel: contains information about transfer function models ')
