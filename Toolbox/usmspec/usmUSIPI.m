function ser = usmUSIPI
%
% series is US IPI
%
data = load(fullfile('data', 'PROJECT_US_MAN_RAW.dat'));
% y = data(:,1);      %1946-I, 2011-III
y = data(57:end, 1); %1960-I, 2011-III
yor = y;
npr = 12; %number of forecasts
freq = 4; % quarterly data
bg_year = 1960;
bg_per = 1;
%outliers detected by TRAMO;    1960-I, 2011-III
%  61 LS    ( 1 1975)
ly = length(y);
Y = zeros(ly+npr, 1); %matrix for regression variables
Y(61:end, 1) = ones(ly+npr-60, 1);
Ycomp = {'level'};


comp.level = [1, 0, 0];
comp.slope = [1, 0.005, NaN];
comp.seas = [1, 0.1, NaN];
comp.irreg = [1, .1, NaN];
comp.cycle = [1, 0.1, NaN];
comp.conout = 'cycle';
twopi = 2 * pi;
comp.cyclep = [0.9, twopi / 20; NaN, NaN];
comp.cycleb = [twopi / 40., twopi / 6.];


ser.yor = yor;
ser.Y = Y;
ser.Ycomp = Ycomp;
ser.bg_year = bg_year;
ser.bg_per = bg_per;
ser.freq = freq;
ser.comp = comp;
ser.olsres = 1;
ser.npr = npr;
ser.lam = -1;
ser.gft = 1;