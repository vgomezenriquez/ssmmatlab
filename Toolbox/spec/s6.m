function ser = s6
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
% THIS IS A SIMULATED SERIES.  THE MODEL IS
%      (1 - B)(1 - B**4)Z(T) = (1 - 0.4B)A(T)
% THERE ARE 178 OBSERVATIONS IN THIS DATA SET.
%

yor = load(fullfile('data', 's6.dat'));
bg_year = 1970;
bg_per = 1;
freq = 4;
Y = [];
ser.yor = yor;
ser.Y = Y;
ser.bg_year = bg_year;
ser.bg_per = bg_per;
ser.freq = freq;
