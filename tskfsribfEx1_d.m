%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example file to illustrate the use of three Kalman filters with a
% multivariate model for 8 U.S. series, 1980-1, 2012-8. The filters are, 1)
% tskf_sribf, 2) tskf plus qr on the OLS model, 3) square root version of
% 2). In this case, the first filter is the fastest because the data are
% multivariate and long.
%
% Model is y_t = X_t*alpha + mu_t + epsilon_t, where mu_t follows the
% model
%
%   mu_{t+1}   = mu_{t}   + K*beta_{t} + eta_{t}
%   beta_{t+1} = beta_{t} +              zeta_{t},
%
%   K = [1  ]      or      K = [ 1   0  ]
%       [b_2]                  [b_2  1  ]
%       [b_3]                  [b_3  c_3]
%        ...                      ...
%       [b_8]                  [b_8  c_8]
%
% Thus, beta_{t} has dimension one or two. Matrix X for the regression
% variables X_{t} is defined in file modstr_viviusa.m, together with the
% other model matrices. There are some missing data.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

%define multivariate model in the following file
tskfsribfExmd

%initial conditions for the filters
ndelta = size(T, 1);
[ins, ii, ferror] = incossm(T, H, ndelta);
chb = 1;

%run three filters. the tskf_sribf is the fastest because the series is
%multivariate and long.
tic
[e1, f1, hd1, Md1, A1, P1, ne1] = tskfsribf(y, X, Z, G, W, T, H, ins, ii, chb);
toc
tic
[e2, f2, hd2, Md2, A2, P2, qyy2, R2, olsres2] = scakfle2(y, X, Z, G, W, T, H, ins, ii, chb);
toc
tic
[e3, f3, hd3, Md3, A3, LP, qyy3, R3, olsres3] = scakflesqrt(y, X, Z, G, W, T, H, ins, ii, chb);
toc
%compare sum of residual sum of squares and determinantal factors
format long g
disp('residual sums of squares: e1''*e1 e2''*e2 e3''*e3')
disp(e1'*e1), disp(e2'*e2), disp(e3'*e3)
disp('determinantal factors: f1 f2 f3')
disp(f1), disp(f2), disp(f3)
format short
