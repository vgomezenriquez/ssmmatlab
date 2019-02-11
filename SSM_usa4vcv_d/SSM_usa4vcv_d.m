%
% Structural model with four variables (USA):
% output, unemployment, inflation and investment.
% The rho parameters are not included in the model.
%
% Reference: ``Estimating Potential Output, Core Inflation
% and the NAIRU as Latent Variables'', by Rafael Domenech
% and Victor Gomez, Journal of Business and Economic Statistics (2006).
%
% Copyright (c) 21 July 2003 by Victor Gomez
% Ministerio de Hacienda, Direccion Gral. de Presupuestos,
% Subdireccion Gral. de Analisis y P.E.,
% Alberto Alcocer 2, 1-P, D-34, 28046, Madrid, SPAIN.
% E-mail: VGomez@sepg.minhap.es
%
% The author assumes no responsibility for errors or damage resulting from the use
% of this code. Usage of this code in applications and/or alterations of it should
% be referenced. This code may be redistributed if nothing has been added or
% removed and no money is charged. Positive or negative feedback would be appreciated.


% usa_2003q1=xlsread('data\usa_2003q1.xls'); % year   lnGDP   dlnp    u   i_y     ifixed_y
usa_2003q1 = load(fullfile('data', 'usa_2003q1.dat'));
% save('data\usa_2003q1.dat','usa_2003q1','-ascii');
%1946:II
s = 4; %seasonality (not used)

[n, nn] = size(usa_2003q1);
m = 4; %number of series
y = zeros(n, m);
y(:, 1) = usa_2003q1(1:n, 2); % Log of GDP
y(:, 2) = usa_2003q1(1:n, 3); % Inflation
y(:, 3) = usa_2003q1(1:n, 4); % Unemployment rate
y(:, 4) = usa_2003q1(1:n, 5); % Gross Private Investment/GDP


figure
subplot(2, 2, 1);
plot(usa_2003q1(1:n, 1), y(:, 1));
title('Log of GDP');

subplot(2, 2, 2);
plot(usa_2003q1(1:n, 1), y(:, 2));
title('Inflation');

subplot(2, 2, 3);
plot(usa_2003q1(1:n, 1), y(:, 3));
title('Unemployment rate');

subplot(2, 2, 4);
plot(usa_2003q1(1:n, 1), y(:, 4));
title('G. P. Investment');

disp('to delete all figures, call function close all')
disp('strike a key to continue')
pause

% initial values for the standard deviations in STAMP:
% level .3679
% slope .2231
% seasonal .1353
% others .6065
%
x0 = [.5, .5, .5, .5, ... %1-4  level disturbances (\sigma_gw, _uw, _xw, _piw)
    .6, .6, .6, ... %5-7  irrg. disturb. (_yw, _uv, _xv), (sigma_piv is concentrated out)
    .5, ... %8     stationary autoregressive parameter: \phi_u
    .5, ... %9     stationary autoregressive parameter: \beta_x
    .5, .1, .1, .1, ... %10-13 stationary autoregressive parameters: \mu_{pi i}
    -.5, .0, .0, ... %14-16 ma parameters (unemployment): \phi_{i}
    .5, -.05, ... %17-18 ma parameters (investment): \beta_{y i}
    .5, ... %19    ma parameters (inflation): \eta_y
    .5, .3142, ... %20-21 cycle parameters:
    .6, ... %22    additional \sigma_yw
    .5, .5, ... %23-24 additional \sigma_piw
    ]
ndelta = 4; %number of diffuse components
nbeta = 4; %number of regression parameters
pfix = [15:16]; %fixed parameters
pvar = [1:14, 17:24]; %free parameters


%parameter optimization


smth = 0;
xv = x0(pvar);
xf = x0(pfix);


mt = 4; %transformation of parameters
if mt > 0
    xv = usatrcv(xv, mt, pfix, pvar, xf);
end

%Levenberg-Marquardt with transformed parameters
info.f = 'usa4vcvf';
info.tr = 1;
info.tolf = 1e-4;
info.tolx = sqrt(info.tolf);
info.max = 300;
info.nu0 = .01;
info.prt = 2;
tic
[x, J, ff, g, iter, conf] = marqdt(info, xv, y, smth, pfix, pvar, xf, 0, mt);
toc


%
% alternative optimization method. Requires the use of the OPTIMIZE toolbox of MATLAB
%
% mt=0; xv=x0(pvar);
% %Options=[];                           %large-scale trust-region method with parameter bounds
% %Options=optimset('Display','iter','Diagnostics','on');
% Options=optimset('Tolfun',1.e-6,'TolX',1.e-6,'Display','iter','Diagnostics','on');
% lb=[]; ub=[];
% tic
%   [x,resnorm,res,exitflag,output,lambda,J]=...
%   lsqnonlin('usa4vcvf',xv,lb,ub,Options,y,s,smth,pfix,pvar,xf,0,mt);
% toc
% exitflag
% output

%end of optimization

if mt == 0
    xut = x;
else
    xut = usautrcv(x, mt, pfix, pvar, xf); %parameters are untransformed
end


xx = x0;
xx(pvar) = xut; %estimated parameters
disp('fixed parameters:'), disp(pfix);


%
% get residuals and estimator
%
smth = -1;
[F, e, hb, Mb, Xt, Pt, nmiss, initf, Pevf] = usa4vcvf(x, y, smth, pfix, pvar, xf, 1, mt);
% save residuals
nname = ['results', filesep, 'residvcv.dat'];
str = ['save ''', nname, ''' e -ascii -double'];
eval(str); % plot residuals
% plot(e)
% pause

% compute residual diagnostics
Ss = e' * e;
Ff = F' * F;
ndrs = (n - 4) * m - nmiss - ndelta; %residual sum of squares
disp('concentrated parameter:')
conp = Ss / (ndrs - nbeta); %estimated sigma square
disp(conp)
disp('square root of concentrated parameter:')
sconp = sqrt(conp);
disp(sconp)
% compute prediction error variance (finite sample)
disp('prediction error variance (finite sample)')
Pevf = Pevf * conp;
disp(Pevf)
disp('standard error (finite sample)')
SPevf = sqrt(diag(Pevf));
disp(SPevf)
%
lagl = round(m*(length(e))^(.5)) + 5;
infr = rescomp(e, lagl, length(pvar), Ss, conp, sconp, Ff, ndrs, nbeta);
%file for output
fname = fullfile('results', 'residvcv.txt');
fid = fopen(fname, 'w');
% fid=1;
%print residual diagnostics
printres(fid, infr);
%close external file
if fid ~= 1
    fclose(fid);
end
%end of residual diagnostics
%outliers
nout = 3;
X = zeros((n - 4)*m, nout);
X(17*m, 1) = 1; %21 = 1951:I   (4 data points are lost)
X(8*m, 2) = 1; %12 = 1948:IV  (4 data points are lost)
for i = 1:30
    X((8 + i)*m, 2) = .7^i; %15 = 1949:III (4 data points are lost)
end
X(19*m, 3) = 1; %23 = 1951:III (4 data points are lost)
% X(11*m,4)=1;  %15 = 1949:III (4 data points are lost)
% for i=1:30
%  X((11+i)*m,4)=.7^i;  %15 = 1949:III (4 data points are lost)
% end
if nout > 0
    Youtg = [];
    t = 1:n - 4;
    for i = 1:nout
        Youtg = [Youtg, X(t, i) * hb(1+i)]; %generate matrix of outlier effects
    end
end
%plot residual diagnostics and outliers
plotres([], [], [], [], [], 1.96, 'residuals', 1, nout, Youtg, 0, [], infr, 1, 1);
disp('to delete all figures, call function close all')
disp('strike a key to continue')
pause

disp('regression parameters')
disp(hb)
disp(Mb)
disp('t-values of regression parameters:')
Mb = Mb * conp;
ttr = hb ./ sqrt(diag(Mb));
disp(ttr)


%
% save results to matlab file param.mat
%
% save param.mat xx xut sconp pfix pvar xf

%standard errors via second derivatives of log likelihood
xutt = xut';
smth = 0;
H = fdhess('logF', xutt, 'usa4vcvf', y, smth, pfix, pvar, xf, 1, 0);
SS = inv(H/2) / ((n - 4) * m - nmiss);
% t-values
disp('t-values:')
tt = zeros(size(xx));
tt(pfix) = NaN;
tt(pvar) = xut ./ sqrt(abs(diag(SS)))';
disp(tt(pvar))

%printing of parameter estimation results

fileprtcv(hb, ttr, xx, tt, sconp)


%filtering
smth = 2;
[F, e, hb, Mb, Xt, Pt, nmiss, initf] = usa4vcvf(x, y, smth, pfix, pvar, xf, 1, mt);
[trendf, strendf, cicf, scicf] = stfilt4(Xt, Pt, xx, y);

%printing of estimated unobserved components
fid = fopen(fullfile('results', 'usa4vcvtf.txt'), 'w');
fprintf(fid, '%15.8f%15.8f%15.8f%15.8f\n', trendf(:, 1:4)');
fclose(fid);
fid = fopen(fullfile('results', 'usa4vcvstf.txt'), 'w');
fprintf(fid, '%15.8f%15.8f%15.8f%15.8f\n', strendf(:, 1:4)');
fclose(fid);

fid = fopen(fullfile('results', 'usa4vcvcf.txt'), 'w');
fprintf(fid, '%15.8f%15.8f%15.8f%15.8f\n', cicf(:, 1:4)');
fclose(fid);
fid = fopen(fullfile('results', 'usa4vcvscf.txt'), 'w');
fprintf(fid, '%15.8f%15.8f%15.8f%15.8f\n', scicf(:, 1:4)');
fclose(fid);


%smoothing
smth = 1;
[F, e, hb, Mb, Xt, Pt, nmiss, initf] = usa4vcvf(x, y, smth, pfix, pvar, xf, 1, mt);
[trends, strends, cics, scics] = stsmt4(Xt, Pt, xx, y, sconp);


%printing of estimated unobserved components

fid = fopen(fullfile('results', 'usa4vcvts.txt'), 'w');
fprintf(fid, '%15.8f%15.8f%15.8f%15.8f\n', trends(:, 1:4)');
fclose(fid);
fid = fopen(fullfile('results', 'usa4vcvsts.txt'), 'w');
fprintf(fid, '%15.8f%15.8f%15.8f%15.8f\n', strends(:, 1:4)');
fclose(fid);

fid = fopen(fullfile('results', 'usa4vcvcs.txt'), 'w');
fprintf(fid, '%15.8f%15.8f%15.8f%15.8f\n', cics(:, 1:4)');
fclose(fid);
fid = fopen(fullfile('results', 'usa4vcvscs.txt'), 'w');
fprintf(fid, '%15.8f%15.8f%15.8f%15.8f\n', scics(:, 1:4)');
fclose(fid);

%plot results. If not desired, uncomment the following line
%return

% confidence interval width: cw=1.69
cw = 1.69;

[nm, m] = size(cicf);
t = 1:nm;

figure
plot(t, cicf(:, 1), t, cicf(:, 2), t, cicf(:, 3), t, cicf(:, 4), ...
    t, zeros(length(t)))
legend('Estimated Cycles: filtering')

figure
plot(t, cicf(:, 1), t, zeros(length(t)), t, cicf(:, 1)+cw*scicf(:, 1), ':', ...
    t, cicf(:, 1)-cw*scicf(:, 1), ':')
legend('Output Gap: filtering')

figure
plot(t, cicf(:, 4), t, zeros(length(t)), t, cicf(:, 4)+cw*scicf(:, 4), ':', ...
    t, cicf(:, 4)-cw*scicf(:, 4), ':')
legend('Cyclical Inflation: filtering')

figure
plot(t, cicf(:, 2), t, zeros(length(t)), t, cicf(:, 2)+cw*scicf(:, 2), ':', ...
    t, cicf(:, 2)-cw*scicf(:, 2), ':')
legend('Cyclical Unemployment: filtering')

figure
plot(t, cicf(:, 3), t, zeros(length(t)), t, cicf(:, 3)+cw*scicf(:, 3), ':', ...
    t, cicf(:, 3)-cw*scicf(:, 3), ':')
legend('Cyclical Investment: filtering')

disp('to delete all figures, call function close all')
disp('strike a key to continue')
pause


figure
plot(t, trendf(:, 1), 'r', t, trendf(:, 1)+cw*strendf(:, 1), ':', ...
    t, trendf(:, 1)-cw*strendf(:, 1), ':', t, y(5:n, 1), 'b')
legend('trend: filtering', '', '', 'Output')

figure
plot(t, trendf(:, 4), 'r', t, trendf(:, 4)+cw*strendf(:, 4), ':', ...
    t, trendf(:, 4)-cw*strendf(:, 4), ':', t, y(5:n, 2), 'b')
legend('trend: filtering', '', '', 'Inflation')

figure
plot(t, trendf(:, 2), 'r', t, trendf(:, 2)+cw*strendf(:, 2), ':', ...
    t, trendf(:, 2)-cw*strendf(:, 2), ':', t, y(5:n, 3), 'b')
legend('trend: filtering', '', '', 'Unemployment')

figure
plot(t, trendf(:, 3), 'r', t, trendf(:, 3)+cw*strendf(:, 3), ':', ...
    t, trendf(:, 3)-cw*strendf(:, 3), ':', t, y(5:n, 4), 'b')
legend('trend: filtering', '', '', 'Investment')

disp('to delete all figures, call function close all')
disp('strike a key to continue')
pause


figure
plot(t, cics(:, 1), t, cics(:, 2), t, cics(:, 3), t, cics(:, 4), ...
    t, zeros(length(t)))
legend('Estimated Cycles: smoothing')

figure
plot(t, cics(:, 1), t, zeros(length(t)), t, cics(:, 1)+cw*scics(:, 1), ':', ...
    t, cics(:, 1)-cw*scics(:, 1), ':')
legend('Output Gap: smoothing')

figure
plot(t, cics(:, 4), t, zeros(length(t)), t, cics(:, 4)+cw*scics(:, 4), ':', ...
    t, cics(:, 4)-cw*scics(:, 4), ':')
legend('Cyclical Inflation: smoothing')

figure
plot(t, cics(:, 2), t, zeros(length(t)), t, cics(:, 2)+cw*scics(:, 2), ':', ...
    t, cics(:, 2)-cw*scics(:, 2), ':')
legend('Cyclical Unemployment: smoothing')

figure
plot(t, cics(:, 3), t, zeros(length(t)), t, cics(:, 3)+cw*scics(:, 3), ':', ...
    t, cics(:, 3)-cw*scics(:, 3), ':')
legend('Cyclical Investment: smoothing')

disp('to delete all figures, call function close all')
disp('strike a key to continue')
pause


figure
plot(t, trends(:, 1), 'r', t, trends(:, 1)+cw*strends(:, 1), ':', ...
    t, trends(:, 1)-cw*strends(:, 1), ':', t, y(5:n, 1), 'b')
legend('trend: smoothing', '', '', 'Output')

figure
plot(t, trends(:, 4), 'r', t, trends(:, 4)+cw*strends(:, 4), ':', ...
    t, trends(:, 4)-cw*strends(:, 4), ':', t, y(5:n, 2), 'b')
legend('trend: smoothing', '', '', 'Inflation')

figure
plot(t, trends(:, 2), 'r', t, trends(:, 2)+cw*strends(:, 2), ':', ...
    t, trends(:, 2)-cw*strends(:, 2), ':', t, y(5:n, 3), 'b')
legend('trend: smoothing', '', '', 'Unemployment')

figure
plot(t, trends(:, 3), 'r', t, trends(:, 3)+cw*strends(:, 3), ':', ...
    t, trends(:, 3)-cw*strends(:, 3), ':', t, y(5:n, 4), 'b')
legend('trend: smoothing', '', '', 'Investment')

disp('to delete all figures, call function close all')
disp('strike a key to continue')
pause


figure
plot(t, cics(:, 1), t, cicf(:, 1), t, zeros(length(t)))
legend('Output gap: sm', 'Output gap: ft')

figure
plot(t, cics(:, 4), t, cicf(:, 4), t, zeros(length(t)))
legend('Cyclical Inflation: sm', 'Cyclical Inflation: ft')

figure
plot(t, cics(:, 2), t, cicf(:, 2), t, zeros(length(t)))
legend('Cyclical Unemployment: sm', 'Cyclical Unemployment: ft')

figure
plot(t, cics(:, 3), t, cicf(:, 3), t, zeros(length(t)))
legend('Cyclical Investment: sm', 'Cyclical Investment: ft')

disp('to delete all figures, call function close all')
disp('strike a key to continue')
pause


figure
plot(t, scics(:, 1), t, scicf(:, 1), t, zeros(length(t)))
legend('ste Output gap: sm', 'ste Output gap: ft')

figure
plot(t, scics(:, 4), t, scicf(:, 4), t, zeros(length(t)))
legend('ste Cyclical Inflation: sm', 'ste Cyclical Inflation: ft')

figure
plot(t, scics(:, 2), t, scicf(:, 2), t, zeros(length(t)))
legend('ste Cyclical Unemployment: sm', 'ste Cyclical Unemployment: ft')

figure
plot(t, scics(:, 3), t, scicf(:, 3), t, zeros(length(t)))
legend('ste Cyclical Investment: sm', 'ste Cyclical Investment: ft')

disp('to delete all figures, call function close all')
disp('strike a key to continue')
pause


figure
plot(t, strends(:, 1), t, strendf(:, 1), t, zeros(length(t)))
legend('ste Output trend: sm', 'ste Output trend: ft')

figure
plot(t, strends(:, 4), t, strendf(:, 4), t, zeros(length(t)))
legend('ste Inflation trend: sm', 'ste Inflation trend: ft')

figure
plot(t, strends(:, 2), t, strendf(:, 2), t, zeros(length(t)))
legend('ste Unemployment trend: sm', 'ste Unemployment trend: ft')

figure
plot(t, strends(:, 3), t, strendf(:, 3), t, zeros(length(t)))
legend('ste Investment trend: sm', 'ste Investment trend: ft')

disp('to delete all figures, call function close all')
