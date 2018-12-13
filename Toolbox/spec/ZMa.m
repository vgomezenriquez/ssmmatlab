function ser = ZMa
%
% Prueba para ver si los residus coinciden con otro orden de regresores
%
ser.autmid = 0;
data = load('data\ZManoB.txt');
yor = data(:, 1);
dr = load('data\ZRManoB.txt');
y1 = dr(:, 1);
y2 = dr(:, 2);
y3 = dr(:, 3);
y4 = dr(:, 4);
y5 = dr(:, 5);
y6 = dr(:, 6);
y7 = dr(:, 7);
y8 = dr(:, 8);
y9 = dr(:, 9);
Y = [y1, y2, y3, y4, y5, y6, y7, y8, y9]; % variables de regresion
rnamesrg = strvcat(['r1'; 'r2'; 'r3'; 'r4'; 'r5'; 'r6'; 'r7'; 'r8'; 'r9']);

bg_year = 1992;
bg_per = 1;
freq = 12;

[ny, my] = size(yor);

ser.lam = 1;
ser.flagm = 0;
ser.p = 3;
ser.dr = 1;
ser.q = 0;
ser.ps = 0;
ser.ds = 1;
ser.qs = 1;
ser.yor = yor;
ser.Y = Y;
ser.rnamesrg = rnamesrg;
ser.bg_year = bg_year;
ser.bg_per = bg_per;
ser.freq = freq;
ser.out = 0;
ser.omet = 1;
ser.C = 3.8;
ser.trad = 0;
ser.easte = 0;
ser.leapy = 0;
ser.gft = 0;
