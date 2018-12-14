function ser = usmbjsgairl
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
% Series is airline series from Box and Jenkins (1976) (series g)


yor = load(fullfile('data', 'bjsgairl.dat'));
bg_year = 1949;
bg_per = 1;
freq = 12;
Y = [];

% ny=size(yor,1);
npr = 24;
% yor=yor(1:ny-npr,:);
comp.level = [1, 0.1, NaN];
comp.slope = [-1, 0., 0];
comp.seas = [1, .1, NaN];
comp.irreg = [1, .1, NaN];


ser.yor = yor;
ser.Y = Y;
ser.bg_year = bg_year;
ser.bg_per = bg_per;
ser.freq = freq;
ser.comp = comp;
ser.npr = npr;
ser.lam = -1;
ser.gft = 1;
