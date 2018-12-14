function ser = bjsaless
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

% Box and Jenkins, (1976)
% SALES DATA WITH LEADING INDICATOR
% THERE ARE 150 OBSERVATIONS
% The identified and estimated model is (SCA Liu, 2005, p. 5.14):
% (1-B)y_t = 0.035 + 4.726*B^3/(1-0.724*B)(1-B)x_t + (1-0.626*B)a_t.

yor = load(fullfile('data', 'bjsales.dat')); %output
n = length(yor); %length of series
x1 = load(fullfile('data', 'bjlead.dat')); %input
bg_year = 1000;
bg_per = 1;
freq = 1;
ninput = 1;
Yin = [];
for i = 1:ninput
    eval(['Yin=[Yin x', num2str(i), '(1:n)];'])
end
yor = yor(1:n);

nlagtf = -1;
tfident = 1;

ser.yor = yor;
ser.Yin = Yin;
ser.bg_year = bg_year;
ser.bg_per = bg_per;
ser.freq = freq;
ser.ninput = ninput;
ser.nlagtf = nlagtf;
ser.tfident = tfident;
