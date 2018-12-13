function ser = turkey_elec
%
% Copyright (c) 21 July 2015 by Victor Gomez
% Ministerio de Hacienda, Direccion Gral. de Presupuestos,
% Subdireccion Gral. de Analisis y P.E.,
% Alberto Alcocer 2, 1-P, D-34, 28046, Madrid, SPAIN.
% Phone : +34-915835439
% E-mail: vgomez@sepg.minhac.es
%
% The author assumes no responsibility for errors or damage resulting from the use
% of this code. Usage of this code in applications and/or alterations of it should
% be referenced. This code may be redistributed if nothing has been added or
% removed and no money is charged. Positive or negative feedback would be appreciated.
%
%Series is Turkey Electricity Demand, data analyzed in De Livera,
%Hyndman and Snyder (2013), Journal of the American Statistical
%Association, ''Forecasting Time Series With Complex Seasonal Patterns
%Using Exponential Smoothing'', 106, 1513-1527.


yor = load(fullfile('data', 'turkey_elec.dat'));
msample = length(yor);
npr = 1096; %number of forecasts
bg_year = 2000;
bg_per = 1;
freq = 7;

%regression variables, two fixed seasonal patterns
modescrs.seas = 2;
modescrs.seasp{1} = [354.37, floor(354.37/2)];
modescrs.seasp{2} = [365.25, floor(365.25/2)];
Y = genfixseaspat(modescrs, msample+npr);

ser.yor = yor;
ser.Y = Y;
ser.lam = 1;
ser.bg_year = bg_year;
ser.bg_per = bg_per;
ser.freq = freq;
ser.npr = npr;
%no automatic model identification, we identify the model
ser.autmid = 0;
ser.dr = 0;
ser.p = 1;
ser.q = 0;
ser.flagm = 1;
%for outlier detection, set ser.out=1
ser.out = 0;
ser.omet = 0;
ser.C = 4.5;
ser.gft = 1;
