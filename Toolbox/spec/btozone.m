function ser = btozone
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
%Series is ozone series from Box and Tiao (1973)


btoz = load(fullfile('data', 'btozone.dat'));
nb = size(btoz, 1);
nyb = nb - 12;
npr = 11;
bg_year = 1955;
bg_per = 1;
freq = 12;
yor = btoz(1:nyb, 1);
Ya = btoz(1:nyb+npr, 2:4);
ct = deltafil(Ya(:, 2:3), 0, 1, 0, 0, freq);
Y = [Ya(:, 1), ct];
rnamesrg = ['int', num2str(1)];
for i = 2:3
    rnamesrg = char(rnamesrg, ['int', num2str(i)]);
end
ser.yor = yor;
ser.Y = Y;
ser.bg_year = bg_year;
ser.bg_per = bg_per;
ser.freq = freq;
ser.rnamesrg = rnamesrg;
ser.npr = npr;
ser.lam = 1;
ser.autmid = 1;
%outlier detection
ser.out = 0;
ser.omet = 1;
ser.C = 3.0;
%graphs
ser.gft = 0;
