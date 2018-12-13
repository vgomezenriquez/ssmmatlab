function ser = tf2
%
% Copyright (c) 21 July 2003 by Victor Gomez
% Ministerio de Hacienda, Direccion Gral. de Presupuestos,
% Subdireccion Gral. de Analisis y P.E.,
% Alberto Alcocer 2, 1-P, D-34, 28046, Madrid, SPAIN.
% Phone : +34-915835439
% E-mail: VGomez@sepg.minhap.es
%
% The author assumes no responsibility for errors or damage resulting from the use
% of this code. Usage of this code in applications and/or alterations of it should
% be referenced. This code may be redistributed if nothing has been added or
% removed and no money is charged. Positive or negative feedback would be appreciated.
%

% simulated trasfer function model (SCA)
% TX2 AND TY2 ARE TWO SIMULATED SERIES.  THE MODEL IS
%      (1 - B)Y(T) = (3.0B - 2.0B**2)(1-B)X(T) + (1 - 0.7B)A(T)
% THERE ARE 130 OBSERVATIONS IN THIS DATA SET.
%input model is
% (1-B)x_t = alpha_t
%
yy = load(fullfile('data', 'vf_tf2.dat'));
x = yy(:, 2); %input
ninput = 1; %number of inputs
yor = yy(:, 1); %output

%fictitious initial date is given
bg_year = 1001;
bg_per = 1;
freq = 1;
nlagtf = -1;
tfident = 1;
Yin = x;

%input forecasts
npr = 12;
ser.npr = npr;
ser.modpred.pred = ones(npr, 1) .* x(end);
%input model (used to compute the mse of the forecasts only)
hr3 = 1;
finv2 = 0;
xd = diferm(x, 1);
[strv, ferror] = estvarmaxpqrPQR(xd, [], freq, [0, 0, 0], [0, 0, 0], hr3, finv2);
ser.modinput.mod = 1;
ser.modinput.phi = [-1., 1.];
ser.modinput.theta = 1.;
ser.modinput.sigma2 = strv.sigmar2;

ser.yor = yor;
ser.Yin = Yin;
ser.bg_year = bg_year;
ser.bg_per = bg_per;
ser.freq = freq;
ser.ninput = ninput;
ser.nlagtf = nlagtf;
ser.tfident = tfident;
ser.gft = 3;
%uncomment the following two lines for no automatic ARIMA model
%identification. The ARIMA model is given by the user
ser.dr = 1;
ser.q = 1;
ser.lam = 1;
ser.pfix = 1;
ser.vfix = -.7;
ser.autmid = 0;
