%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example file to illustrate the use of three Kalman filters with a
% univariate model for the airline series from Box and Jenkins (1976). The
% filters are, 1) tskf_sribf, 2) tskf plus qr on the OLS model, 3) square
% root version of 2). In this case, the second filter is the fastest
% because the data are univariate and short.
%
% Model is univariate structural model with trend, slope, trigonometric
% seasonality, cycle, irregular and autoregressive component. The slope is
% fixed and there are some missing data. Data are monthly.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

%define and estimate univariate model in the following file
tskfsribfExud

%initial conditions for the filters
ndelta = size(T, 1);
[ins, ii, ferror] = incossm(T, H, ndelta);
chb = 1;

%run three filters. the tskf_sribf is the fastest because the series is
%multivariate and long.
tic
[e1, f1, hd1, Md1, A1, P1, ne1] = tskfsribf(yl, X, Z, G, W, T, H, ins, ii, chb);
toc
tic
[e2, f2, hd2, Md2, A2, P2, qyy2, R2, olsres2] = scakfle2(yl, X, Z, G, W, T, H, ins, ii, chb);
toc
tic
[e3, f3, hd3, Md3, A3, LP, qyy3, R3, olsres3] = scakflesqrt(yl, X, Z, G, W, T, H, ins, ii, chb);
toc
%compare sum of residual sum of squares and determinantal factors
format long g
disp('residual sums of squares: e1''*e1 e2''*e2 e3''*e3')
disp(e1'*e1), disp(e2'*e2), disp(e3'*e3)
disp('determinantal factors: f1 f2 f3')
disp(f1), disp(f2), disp(f3)
format short
