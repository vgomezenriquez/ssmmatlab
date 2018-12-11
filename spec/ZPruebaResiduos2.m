function ser = ZPruebaResiduos2
%
% Prueba para ver si los residus coinciden con los de TSW y JDemetra+
%
ser.autmid = 0;
data = load('data\ZPruebaResiduos.txt');
yor = data(:, 1);
y1 = data(:, 2); %number of kilometers driven
Y = [y1]; %matrix for regression variables
rnamesrg = 'regresor';

bg_year = 2006;
bg_per = 1;
freq = 4;

[ny, my] = size(yor);


ser.lam = 1;
ser.dr = 0;
ser.ds = 1;
ser.p = 2;
ser.ps = 0;
ser.q = 0;
ser.qs = 0;
ser.flagm = 0;
ser.yor = yor;
ser.Y = Y;
ser.rnamesrg = rnamesrg;
ser.bg_year = bg_year;
ser.bg_per = bg_per;
ser.freq = freq;
ser.out = 0;
ser.omet = 1;
ser.C = 3.8;
ser.gft = 1;
