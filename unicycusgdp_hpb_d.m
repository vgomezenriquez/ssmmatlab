%********************************************************************
%       UNIVARIATE CYCLE ESTIMATION OF US GDP
%                1953.Q2 - 2007.Q3
%
%  METHOD: HODRICK-PRESCOTT FILTER AND BAND-PASS FILTER
%
%  Usually, the HP filter is applied without a model for the series. When
%  this is attempted with the band-pass filter, it does not work. But, if
%  this last filter is applied to the series detrended with the HP filter,
%  it works.
%
%*********************************************************************

clear;

% application of the HP filter.
out = usmestos('usmUSIPIHP');
y = out.orig;
s = out.freq; %seasonal frequency
bg_year = out.nziyip(2); %initial year
bg_per = out.nziyip(3); %initial period
datei = cal(bg_year, bg_per, s);

% Read the Excel file with the NBER recession dates
% recdates = xlsread('data\NBER_chron_short.xls');
recdates = load(fullfile('data', 'NBER_chron_short.dat'));
[mrec, nrec] = size(recdates);

% HP trend
thp = out.model.StochCc(:, 1);
% HP Cycle
chp = y - thp;

% Plot the HP cycle with the NBER recession dates
% Setting fot the plots
ydigit = 'yyyy';
nberrecplot(recdates, chp)
hold on
tsplot(chp, datei, 'GDP HP Cycle', ydigit)
hold off
pause

%************************************************************************
%                    BAND-PASS FILTER
%
% Design a band-pass filter to obtain a well defined cycle and a
% relatively smooth trend. Frequencies are expressed divided by pi.

D(1) = .1;
D(2) = .1;
xp1 = .0625;
xp2 = .3;
xs = .4;
% Tangent band-pass filter
[compbp, ferror] = dbptanbut(D, xp1, xp2, xs);
% plot gain function of the tangent band-pass filter
figure
ggbptanbut(D, xp1, xp2, xs, compbp.Di, compbp.Alph, compbp.Lambda)
pause

%Application of the band-pass filter without a model to the GDP series. To
%this end, we use the model implied by the band-pass filter:
%
%  y_t = s_t + n_t,
%
%where s_t follows the model
%
% (1 - 2*\alpha*B + B^2)^d s_t = (1 - B^2)^d b_t,
%
%\sigma^2_b = 1, and \sigma^2_n = Lambda. Here, d=compbp.Di, \alpha =
%compbd.Alph and Lambda = compbp.Lambda.
%We use the cascade form for the filter for numerical reasons.
%
alpha = compbp.Alph;
den = [1., -2 * alpha, 1.];
Alpha = [-1., 0, 1.];
Jp = 1.;
[Tp, Hp, Zp, ferror1] = akaikessm1(den, Alpha);
[Tsp, Hsp, Zsp, Jsp, ferror] = cascadessm1(Tp, Hp, Zp, Tp, Hp, Zp, Jp);
for i = 1:compbp.Di - 3
    [Tsp, Hsp, Zsp, Jsp, ferror] = cascadessm1(Tsp, Hsp, Zsp, Tp, Hp, Zp, Jp);
end
[Tp, Hp, Zp, ferror1] = akaikessm2(den, Alpha);
Jp = 0;
[Tsp, Hsp, Zsp, Jsp, ferror] = cascadessm1(Tsp, Hsp, Zsp, Tp, Hp, Zp, Jp);
X = [];
W = [];
T = Tsp;
Z = Zsp;
H = [Hsp, zeros(size(T, 1), 1)];
G = [0, sqrt(compbp.Lambda)];
%number of unit roots in Tsp. Equal to 2 times Di.
ndelta = 2 * compbp.Di;
[ins, ii, ferror] = incossm(T, H, ndelta);
[KKP, PT, a, b] = scakfssqrt(y, X, Z, G, W, T, H, ins, ii);
plot(KKP(:, 1))
legend('filter model: a mess')
disp('strike any key to continue')
pause
%A mess because the model implied by the filter is very different from the
%model followed by the IPI series. Therefore, it is not a good idea to use
%a band-pass filter without a model. At least, not a band-pass filter with
%a model-based interpretation.

%But if we apply the band-pass filter to the "detrended series" (using the
%HP filter), now it works.
[KKP, PT, a, b] = scakfssqrt(chp, X, Z, G, W, T, H, ins, ii);
tp = 1:size(y, 1);
cyclebp = KKP(:, 1);
plot(tp, chp, tp, cyclebp, 'r');
legend('cycle hp', 'band-pass')
disp('strike any key to continue')
pause
close all
