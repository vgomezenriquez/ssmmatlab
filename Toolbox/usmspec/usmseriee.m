function ser = usmseriee
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
% Series is series e from Box and Jenkins (1976)
%  WOLFER SUNSPOT NUMBERS: YEARLY, 1770 - 1869


yor = load(fullfile('data', 'seriee.dat'));
bg_year = 1770;
bg_per = 1;
freq = 1;
Y = [];

npr = 10;
%define univariate structural model: trend, slope, trigonometric
%seasonality, cycle, irregular and autoregressive component
comp.level = [1, 0.1, NaN];
% comp.slope=[-1 0. 0];
% comp.irreg=[1 .1  NaN];
comp.cycle = [1, .1, NaN];
twopi = 2 * pi;
comp.cyclep = [.9, twopi / 10.; NaN, NaN];
comp.cycleb = [twopi / 15., twopi / 5.];
comp.freq = freq;
% %create some missing values in the series
% yor(2:7)=NaN(6,1);
% yor(20)=NaN;


ser.yor = yor;
ser.Y = Y;
ser.bg_year = bg_year;
ser.bg_per = bg_per;
ser.freq = freq;
ser.comp = comp;
ser.npr = npr;
ser.lam = -1;
ser.gft = 1;
