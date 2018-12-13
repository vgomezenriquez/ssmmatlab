function ser = Seatbelt_arima
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

data = load(fullfile('data', 'Seatbelt.dat'));
yor = data(:, 1);
y1 = data(:, 4); %number of kilometers driven
y2 = data(:, 5); %price of oil
Y = [y1, y2]; %matrix for regression variables
rnamesrg = strvcat(['km driven'; ... %names for regression variables
    'oil      ']);

bg_year = 1969;
bg_per = 1;
freq = 12;

ser.yor = yor;
ser.Y = Y;
ser.rnamesrg = rnamesrg;
ser.bg_year = bg_year;
ser.bg_per = bg_per;
ser.freq = freq;
ser.out = 1;
ser.omet = 1;
ser.C = 3.0;
ser.gft = 1;
