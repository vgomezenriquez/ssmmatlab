%Example of automatic model identification and estimation of a regression
%model with ARIMA  errors. The program performs automatic outlier
%detection.
%Series is car drivers killed or seriously injured in Great Britain from
%January 1969 to December 1984 (Durbin and Koopman, 2012).
%Two explanatory variables are included in the model, the price of oil and
%the number of kilometers driven.
%
%The specification file, Seatbelt_arima.m, is in the subdirectory spec. In
%this specification file, the instructions for the program are given. The
%default values for the program are in the script file arimadefval.m.
%

out = arimaestos('Seatbelt_arima');

disp(' ')
disp('Details on the identified and estimated model are in the')
disp('file "Seatbelt_arima.txt" in the subdirectory "results"')
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