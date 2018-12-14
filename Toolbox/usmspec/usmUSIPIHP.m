function ser = usmUSIPIHP
%
% series is US IPI. The structural model is the one corresponding to the
% Hodrick-Prescott filter for quarterly series (lambda = 1600)
%
% Read the data
data = load(fullfile('data', 'FINAL_CKZ_JAE2.dat'));
yy = data(:, 2:end);
%calendar for the whole time span 1947.M1 - 2007.M9
bgt_year = 1947;
bgt_per = 1;
freqm = 12;
dateit = cal(bgt_year, bgt_per, freqm);
%initial date for the data
idate = ical(1953, 6, dateit);
%final date for the data
fdate = ical(2007, 9, dateit);

yy = yy(idate:fdate, :);
y = yy(:, 6); % GDP
% Remove NaNs
y(isnan(y)) = [];

%calendar for the considered time span 1953.Q2 - 2007.Q3
bg_year = 1953;
bg_per = 2;
freqq = 4;

Y = []; %matrix for regression variables
%-----------------------------------------------------------------
%                   CYCLE ESTIMATION (Hodrick-Prescott)
%

% Specify components, all parameters are fixed. Note that the irregular
% seems to be not fixed. But it will be concentrated out and, therefore, it
% will be made fixed.
comp.level = [1, 0, 0];
lambda = 1600;
sdsl = sqrt(1/lambda); % Restriction for the slope standard deviation
comp.slope = [1, sdsl, 0];
comp.irreg = [1, .1, NaN];
comp.freq = freqq;
npr = 0;

ser.yor = y;
ser.Y = Y;
ser.bg_year = bg_year;
ser.bg_per = bg_per;
ser.freq = freqq;
ser.comp = comp;
ser.olsres = 1;
ser.npr = npr;
ser.lam = 1;
ser.nlestim = 0;
ser.gft = 1;