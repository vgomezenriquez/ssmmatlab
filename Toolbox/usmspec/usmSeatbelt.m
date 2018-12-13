function ser = usmSeatbelt
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
%Example of estimation of a univariate structural model
%Series is car drivers killed or seriously injured in Great Britain from
%January 1969 to December 1984 (Durbin and Koopman, 2012).
%Two explanatory variables are included in the model, the price of oil and
%the number of kilometers driven.
data = load(fullfile('data', 'Seatbelt.dat'));
y = data(:, 1);
y1 = data(:, 4); %number of kilometers driven
y2 = data(:, 5); %price of oil
yor = y;
lam = 1;
Y = [y1, y2]; %matrix for regression variables
Ycomp = {'level', 'level'};
npr = 0; %number of forecasts
%define univariate structural model: trend, trigonometric
%seasonality, and irregular component
comp.level = [1, 0.1, NaN];
comp.seas = [2, .1, NaN];
comp.irreg = [1, .1, NaN];
freq = 12;
bg_year = 1969;
bg_per = 1;

ser.yor = yor;
ser.Y = Y;
ser.Ycomp = Ycomp;
ser.bg_year = bg_year;
ser.bg_per = bg_per;
ser.freq = freq;
ser.comp = comp;
ser.npr = npr;
ser.lam = lam;
ser.gft = 1;
