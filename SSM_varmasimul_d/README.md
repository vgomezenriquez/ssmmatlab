[<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/banner.png" width="888" alt="Visit QuantNet">](http://quantlet.de/)

## [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/qloqo.png" alt="Visit QuantNet">](http://quantlet.de/) **SSM_varmasimul_d** [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/QN2.png" width="60" alt="Visit QuantNet 2.0">](http://quantlet.de/)

```yaml

Name of QuantLet: SSM_varmasimul_d

Published in: Linear Time Series With MATLAB and Octave

Description: 'Firstly, a time series following a VARMA model is simulated. Then, the model is estimated by the Hannan-Rissanen method.'

Keywords: time-series, VARMA model, simulation, estimation, exact maximum likelihood method

Author: Víctor Gómez

Submitted: Wed, December 19 2018 by Víctor Gómez


```

### MATLAB Code
```matlab

%script file to simulate a VARMA model
%
% Let the VARMA model be
%
%   Phi(B) y_t = Theta(B) a_t,
%
% where B is the backshift operator, By_t = y_{t-1}.
%

clear
l = 50; %number of initial observations to be discarded
N = 150; %number of observations of the simulated series
s = 2; %number of variables
seed = 20; %this is to always generate the same series

%polynomial matrices Phi and Theta
Phi(:, :, 1) = eye(s);
Phi(:, :, 2) = -[0.7, 0.0; 0.4, 0.9];
Theta(:, :, 1) = eye(s);
Theta(:, :, 2) = -[0.5, 0.4; 0.0, 0.5];

%covariance matrix of the a_t innovartions
S = [1.0, 1.2; 1.2, 4.0];

%simulate Y_t
[y, ferror] = varmasim(l, N, Phi, Theta, S, seed);

%estimate model using HR method
seas = 1;
x = [];
[strv, ferror] = estvarmaxpqrPQR(y, x, seas, [1, 1, 0], [0, 0, 0]);

disp(' ');
disp('***** Estimated Model  *****');
disp(' ');
clear in
in.fid = 1;
tit = 'Phi';
strt = 1;
mprintar(strv.phis(:, :, 2), in, tit, strt);
disp(' ')
tit = 'Theta';
strt = 1;
mprintar(strv.thetas(:, :, 2), in, tit, strt);
disp(' ')
tit = 'Sigma';
mprintar(strv.sigmar2, in, tit);
disp(' ')
disp(' ')
disp('***** T-values  *****');
disp(' ');
clear in
in.fid = 1;
tit = 'tv-Phi';
strt = 1;
mprintar(strv.phitv(:, :, 2), in, tit, strt);
disp(' ')
tit = 'tv-Theta';
strt = 1;
mprintar(strv.thetatv(:, :, 2), in, tit, strt);

%setup model
Phis = eye(s);
Ths = eye(s);
Phi = strv.phis;
Theta = strv.thetas;
Sigma = strv.sigmar2;
freq = 1;

%create structure and put model into state space form
[str, ferror] = suvarmapqPQ(Phi, Theta, Phis, Ths, Sigma, freq);
%matrix for regression varialbes
Y = [];

%estimate model
result = varmapqPQestim(y, str, Y);

%estimated and fixed parameters
xvf = result.xvf;
xf = result.xf;

%create estimated model
[phif, thf, Phif, Thf, Lf, ferror] = pr2varmapqPQ(xvf, xf, str);

%t-values
tvf = result.tv;

%create matrices with t-values
[phitvf, thtvf, Phitvf, Thtvf, Ltvf, ferror] = pr2varmapqPQ(tvf, xf, str);

disp(' ');
disp('***** Estimated Model  *****');
disp(' ');
clear in
in.fid = 1;
tit = 'Phi';
strt = 1;
mprintar(phif(:, :, 2), in, tit, strt);
disp(' ')
tit = 'Theta';
strt = 1;
mprintar(thf(:, :, 2), in, tit, strt);
disp(' ')
tit = 'Sigma';
mprintar(result.Sigmar, in, tit);
disp(' ')

disp('***** t-values  *****');
disp(' ');
clear in
in.fid = 1;
tit = 'tv-Phi';
strt = 1;
mprintar(phitvf(:, :, 2), in, tit, strt);
disp(' ')
tit = 'tv-Theta';
strt = 1;
mprintar(thtvf(:, :, 2), in, tit, strt);
disp(' ')
tit = 'tv-L';
mprintar(Ltvf, in, tit);

```

automatically created on 2019-02-11