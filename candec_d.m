%script file to perform a canonical decomposition of an ARIMA model

clear
%graphs for spectra
fplot = 1;

%model
s = 12;
dr = 1;
ds = 1;

%AR
phi(:, :, 1) = 1.; %regular part
% phi(:,:,2)=-.95;
Phi(:, :, 1) = 1.; %seasonal part
Phi(:, :, 2) = -.9;

%MA
th(:, :, 1) = 1; %regular part
th(:, :, 2) = -.4;
% th2(:,:,1)=1.; th2(:,:,2)=-.3; th=pmatmul(th,th2);
Th(:, :, 1) = 1.; %seasonal part
Th(:, :, 2) = -.6;
%standard deviation of the innovations
sconp = 1.;

% set up trend-cycle and seasonal polynomials for the canonical
% decomposition
[phir, phis, thr, ths, phirst] = arima2rspol(phi, Phi, th, Th, s, dr, ds);


%perform canonical decomposition
[comp, ierrcandec] = candec(phir, phis, thr, ths, phirst, s, dr, ds, sconp);

sigma2i = comp.itvar; % variance of the irregular component (*)
%(*) in units of the series model innovations

if sigma2i < 0
    disp('irregular spectrum negative, sigma2i=')
    disp(sigma2i)
    pause
end

%plot spectra
if fplot == 1
    out.compcd = comp;
    plotspcd(out)
end

thrc = comp.ptnum; % trend-cycle numerator
phir = comp.ptden; % trend-cycle denominator
%   comp.ptnur;    % number of nonstationary roots in phir
sigma2r = comp.ptvar; % variance of the trend-cycle innovations (*)
thsc = comp.stnum; % seasonal numerator
phis = comp.stden; % seasonal denominator
%   comp.stnur; % number of nonstationary roots in phis
sigma2s = comp.stvar; % variance of the seasonal innovations (*)
thtc = comp.rt; % transitory component (MA term)
sigma2t = comp.rtvar; % variance of the transitory component innovations (*)
sigma2i = comp.itvar; % variance of the irregular component (*)
phitc = comp.phi; % stationary AR trend polynomial
%(*) in units of the series model innovations

disp('trend-cycle numerator:')
disp(thrc)
disp('trend-cycle denominator:')
disp(phir)
disp('variance of the trend-cycle innovations (*)')
disp(sigma2r)
pause
disp('seasonal numerator:')
disp(thsc)
disp('seasonal denominator:')
disp(phis)
disp('variance of the seasonal innovations (*)')
disp(sigma2s)
pause
disp('stationary AR trend polynomial')
disp(phitc)
pause
disp('variance of the irregular component (*)')
disp(sigma2i)
if ~isempty(thtc)
    disp('transitory numerator:')
    disp(thtc)
    phit = 1.;
    disp('transitory denominator:')
    disp(phit)
    disp('variance of the transitory innovations (*)')
    disp(sigma2t)
end
disp('(*) in units of var(A)')
