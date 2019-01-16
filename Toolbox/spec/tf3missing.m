function ser = tf3missing
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
% TX3 AND TY3 ARE TWO SIMULATED SERIES.  THE MODEL IS
%      (1 - B)Y(T) = (4.0B**2 - 4.0B**3 + 1.0B**4)(1 - B)X(T)
%                    + (1 - 0.7B)A(T)
% THERE ARE 115 OBSERVATIONS IN THIS DATA SET.
nombre = load(fullfile('data', 'vf_tf3.dat'));
x = nombre(:, 2); %input
ninput = 1; %number of inputs
yor = nombre(:, 1); %output

%fictitious initial date is given
bg_year = 1001;
bg_per = 1;
freq = 1;
nlagtf = -1;
tfident = 1;
Yin = x;

%Here, we specify some missing values
yor(2) = NaN;
yor(34:40) = NaN(7, 1);
yor(end-3) = NaN;

%input fictitious forecasts
npr = 12;
ser.npr = npr;
twxn = 2 * x(end);
ser.modpred(1).pred = zeros(npr, 1);
ser.modinput.mod = 0;
for i = 1:npr
    ser.modpred(1).pred(i) = twxn - x(end-i); %use symmetry to obtain the forecasts of x1
end

ser.yor = yor;
ser.Yin = Yin;
ser.bg_year = bg_year;
ser.bg_per = bg_per;
ser.freq = freq;
ser.ninput = ninput;
ser.nlagtf = nlagtf;
ser.tfident = tfident;
ser.gft = 1;
ser.olsres = 1;
