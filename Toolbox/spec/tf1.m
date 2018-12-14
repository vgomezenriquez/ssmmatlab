function ser = tf1
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
% TX1 AND TY1 ARE TWO SIMULATED SERIES.  THE MODEL IS
%      (Y(T) - YBAR) = (3.0B**2)/(1 - 0.5B)(X(T) - XBAR) + A(T)
% WHERE XBAR AND YBAR ARE THE MEANS OF THE X(T) AND Y(T) SERIES
% THERE ARE 125 OBSERVATIONS IN THIS DATA SET.
nombre = load(fullfile('data', 'vf_tf1.dat'));
x = nombre(:, 2); %input
ninput = 1; %number of inputs
yor = nombre(:, 1); %output

%subtract mean from input and output
x = x - mean(x);
yor = yor - mean(yor);

%fictitious initial date is given
bg_year = 1001;
bg_per = 1;
freq = 1;
nlagtf = -1;
tfident = 1;
Yin = x;

%input fictitious forecasts
npr = 12;
ser.npr = npr;
ser.modinput.mod = 0;
twxn = 2 * x(end);
ser.modpred.pred = zeros(npr, 1);
for i = 1:npr
    ser.modpred.pred(i) = twxn - x(end-i); %use symmetry to obtain the forecasts of x1
end

ser.yor = yor;
ser.Yin = Yin;
ser.bg_year = bg_year;
ser.bg_per = bg_per;
ser.freq = freq;
ser.ninput = ninput;
ser.lam = -1;
ser.nlagtf = nlagtf;
ser.tfident = tfident;
ser.gft = 1;
%uncomment the following two lines for no automatic ARIMA model
%identification. The ARIMA model is given by the user (white noise)
% ser.dr=0;
% ser.pfix=1; ser.vfix=0.; ser.autmid=0;
%uncomment the following lines for no automatic TF model identification.
%the TF model is given by the user
% ser.tfident=0;
% ser.delay=2; ser.ma=0; ser.ar=1;
