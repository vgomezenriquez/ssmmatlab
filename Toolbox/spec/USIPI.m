function ser = USIPI
%
% series is US IPI
%
data = load(fullfile('data', 'PROJECT_US_MAN_RAW.dat'));
% y = data(:,1);      %1946-I, 2011-III
y = data(57:end, 1); %1960-I, 2011-III
yor = y;
freq = 4; % quarterly data
bg_year = 1960;
bg_per = 1;
Y = []; %matrix for regression variables
npr = 0; %number of forecasts


ser.lam = 1;
ser.dr = 1;
ser.ds = 1;
ser.p = 1;
ser.ps = 0;
ser.q = 0;
ser.qs = 1;
ser.flagm = 0;
ser.yor = yor;
ser.Y = Y;
ser.npr = npr;
ser.autmid = 0;

ser.bg_year = bg_year;
ser.bg_per = bg_per;
ser.freq = freq;
ser.out = 1;
ser.omet = 1;
ser.C = 4.5;
ser.gft = 1;
