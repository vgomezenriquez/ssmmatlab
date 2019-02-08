%Example of automatic model identification and estimation of several ARIMA
%models and one transfer function model.
%The specification file, metafile.txt, is in the subdirectory spec. In this
%specification file, the list with the individual models is given. The
%default values for the program are in the script file arimadefval.m.
%

fmeta = 1;
out = arimaestos('metafile', fmeta);

disp(' ')
disp('Details on the identified and estimated models are in the')
disp('corresponding files ".txt" in the subdirectory "results"')
disp(' ')
disp('The file "summary.txt" in the directory "results" contains')
disp('a summary of the results')
