function ser = usmcgwage
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
% Series is German GDP real wage, seasonally adjusted

data = load(fullfile('data', 'PROJECTDATA.dat'));
data(any(isnan(data)'), :) = [];

yor = data(:, 4);
Y = []; %matrix for regression variables
npr = 10; %number of forecasts
freq = 4; % quarterly data
bg_year = 1970;
bg_per = 1;

% Specify components, initial values and fixed parameters
comp.level = [1, 0, 0];
comp.slope = [1, .005, NaN];
comp.irreg = [1, 0.1, NaN];
comp.cycle = [1, .1, NaN];
twopi = 2 * pi;
comp.cyclep = [0.9, twopi / 40; NaN, NaN];
comp.cycleb = [twopi / 60., twopi / 6.];

ser.yor = yor;
ser.Y = Y;
ser.bg_year = bg_year;
ser.bg_per = bg_per;
ser.freq = freq;
ser.comp = comp;
ser.npr = npr;
ser.lam = 1;
ser.gft = 1;
