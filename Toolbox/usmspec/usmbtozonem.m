function ser = usmbtozonem
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
%Example of estimation of a univariate structural model with missing values
%Series is ozone series from Box and Tiao (1973)
%

btoz = load(fullfile('data', 'btozone.dat'));
nb = size(btoz, 1);
nyb = nb - 12;
yor = btoz(1:nyb, 1);
bg_year = 1955;
bg_per = 1;
freq = 12;
Ya = btoz(:, 2:4);
ct = deltafil(Ya(:, 2:3), 0, 1, 0, 0, freq);
Ya = [Ya(:, 1), ct];
Y = Ya(1:nb, :);
npr = 12; %number of forecasts

%create some missing values in the series
yor(2:7) = NaN(6, 1);
yor(20) = NaN;


%define univariate structural model: trend, slope, trigonometric
%seasonality, cycle, irregular and autoregressive component
comp.level = [-1, 0, 0];
% comp.slope=[-1 0. 0];
comp.seas = [2, .1, NaN];
% comp.irreg=[1 .1  NaN];
comp.ar = [1, .1, NaN];
comp.arp = [-.1; NaN];
comp.freq = freq;

ser.yor = yor;
ser.Y = Y;
ser.bg_year = bg_year;
ser.bg_per = bg_per;
ser.freq = freq;
ser.comp = comp;
ser.npr = npr;
ser.lam = -1;
ser.gft = 1;
