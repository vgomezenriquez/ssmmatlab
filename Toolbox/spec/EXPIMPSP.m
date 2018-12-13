function ser = EXPIMPSP
%
% Example series from TSW+
%
y = load(fullfile('data', 'EXPIMPSP.dat'));
yor = y;
freq = 12; % monthly data
bg_year = 1976;
bg_per = 1;
Y = [];


ser.lam = -1;
% ser.dr=0; ser.ds=0;
% ser.p=2; ser.ps=1;
% ser.q=1; ser.qs=0;
% ser.flagm=1;
ser.yor = yor;
ser.Y = Y;

ser.bg_year = bg_year;
ser.bg_per = bg_per;
ser.freq = freq;
ser.autmid = 1;
ser.out = 1;
ser.omet = 1;
ser.C = 3.;
ser.gft = 1;
