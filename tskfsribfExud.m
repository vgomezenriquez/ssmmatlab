% load data.
y = load(fullfile('data', 'bjsgairl.dat'));
x = [];
yl = log(y); %transform series
Y = []; %matrix for regression variables

%define univariate structural model: trend, slope, trigonometric
%seasonality, cycle, irregular and autoregressive component
comp.level = [1, 0.1, NaN];
comp.slope = [-1, 0, 0];
comp.seas = [2, .1, NaN];
comp.irreg = [1, .1, NaN];
freq = 12;
comp.freq = freq;
bg_year = 1949;
bg_per = 1;
datei = cal(bg_year, bg_per, freq);
comp.datei = datei;
npr = 0;

%create structure and put model into state space form
[str, ferror] = suusm(comp, yl, Y, npr);

%estimate model
[result, str] = usmestim(yl, str);

xv = result.xvf;
xf = result.xf;
sigmac = sqrt(result.sigma2c);
SPevf = result.SPevf;
[X, Z, G, W, T, H, ins, ii, ferror] = pr2usm(xv, xf, str);
%Pevf (innovations variance) is concentrated out
G = (G * sigmac) / SPevf;
H = (H * sigmac) / SPevf;
