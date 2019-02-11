[<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/banner.png" width="888" alt="Visit QuantNet">](http://quantlet.de/)

## [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/qloqo.png" alt="Visit QuantNet">](http://quantlet.de/) **SSM_arimatf9_d** [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/QN2.png" width="60" alt="Visit QuantNet 2.0">](http://quantlet.de/)

```yaml

Name of QuantLet: SSM_arimatf9_d

Published in: Linear Time Series With MATLAB and Octave

Description: 'Several ARIMA and transfer function models are handled in one single run of the program.'

Keywords: time-series, ARIMA model, transfer function model, estimation, automatic model identification

Author: Víctor Gómez

Submitted: Tue, January 10 2019 by Víctor Gómez

```

### MATLAB Code
```matlab

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

```

automatically created on 2019-02-11