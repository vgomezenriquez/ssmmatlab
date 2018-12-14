%
% Example of ''Bivariate structural time series analysis'' in the book by
% Durbin and Koopman (2201), p. 167. There is a discrepancy with the same
% example in the second edition of this book (2012), p. 195. I think the
% correct results are those of the first edition.
%
%
% Copyright (c) 21 July 2003 by Victor Gomez
% Ministerio de Hacienda, Direccion Gral. de Presupuestos,
% Subdireccion Gral. de Analisis y P.E.,
% Alberto Alcocer 2, 1-P, D-34, 28046, Madrid, SPAIN.
% Phone : +34-915835439
% E-mail: VGomez@sepg.minhap.es
%
% The author assumes no responsibility for errors or damage resulting from the use
% of this code. Usage of this code in applications and/or alterations of it should
% be referenced. This code may be redistributed if nothing has been added or
% removed and no money is charged. Positive or negative feedback would be appreciated.
%

clear
data = load(fullfile('data', 'Seatbelt.dat'));


nb = 24;
y1 = data(1:end-nb, 2);
y2 = data(1:end-nb, 3);
y = [y1, y2];

%regression variables: km driven and price of oil
Y = data(1:end-nb, 4:5);

%plot both series
bg_yearm = 1969;
bg_perm = 1;
freqm = 12;
dateim = cal(bg_yearm, bg_perm, freqm);
tsplot(y1, dateim, 'front seats');
disp('strike a key to continue')
pause
tsplot(y2, dateim, 'rear seats');
disp('strike a key to continue')
pause
close all


vnames = char('front seats ','rear seats');
tsplot(y, dateim, vnames);
disp('strike a key to continue');
pause
close all

tname = 'frontrear';
fname = fullfile('results', 'frontrear.txt');

%first, form univariate structural model and, from this, construct
%bivariate structural model (this is easier than constructing the bivariate
%structural model from scratch)

%univariate structural model defined according to ssm_matlab
YY = [];
comp.level = [1, .1, NaN];
comp.seas = [2, .1, NaN];
comp.irreg = [1, .1, NaN];
comp.conout = 'level';
comp.freq = freqm;
comp.datei = dateim;
npr = 0;
[models, ferror] = suusm(comp, y1, YY, npr);
%end of univariate structural model defined according to ssm_matlab

nr = length(models.pvar);

% Bivariate structural model
nvar = 2;
Tb = kron(models.T, eye(nvar));
Hb = kron(models.H, eye(nvar));
Zb = kron(models.Z, eye(nvar));
Gb = kron(models.G, eye(nvar));
insb = kron(models.ins, eye(nvar));
[nt, mt] = size(models.T);
ib = [0, 0, 0, nvar * nt];
Yb = kron(Y, eye(nvar));

nbeta = size(Yb, 2);

% order in the standard deviations (conc):
% 1 2 3    - level
% 4 5 6    - slope
% 7 8 9    - seasonal
% 10 11 12 - irregular
x0 = [.1, .1, ... %sigma level s21 s22 (s11 one parameter is concentrated out)
    .1, .1, .1, ... %sigma seasonal s11 s21 s22;
    .1, .1, .1, ... %sigma irreg. s11 s21 s22
    ];


% pfix=[ ];             %fixed parameters
% pvar=[1:8 ];       %free parameters
pfix = []; %fixed parameters
pvar = [1:8]; %free parameters


xv = x0(pvar);
xf = x0(pfix);

stordb = [2, 3, 4, 5, 6, 7, 8, 9]; %standard deviations to be estimated
concb = 1; %standard deviation concentrated out
[nhb, mhb] = size(Hb);
xp = zeros(1, 9);
xp(stordb) = x0;
xp(concb) = 1; %concentrated parameter
Hb(1, 1) = xp(1);
Hb(2, 1) = xp(2);
Hb(2, 2) = xp(3);
j = 5;
for i = 4:2:nhb
    Hb(i, j) = xp(4);
    Hb(i-1, j) = xp(5);
    Hb(i, j+1) = xp(6);
    j = j + 2;
end
Gb(1, mhb-1) = xp(7);
Gb(2, mhb-1) = xp(8);
Gb(2, mhb) = xp(9);

models.X = Yb;
models.T = Tb;
models.H = Hb;
models.Z = Zb;
models.G = Gb;
models.i = ib;
models.ins = insb;
models.stord = stordb;
models.conc = concb;
models.pvar = pvar;
models.pfix = pfix;
models.x = x0;
models.xv = xv;
models.xf = xf;
% end of bivariate structural model

clear Tb Hb Zb Gb ib insb stordb concb;


s = freqm;
chb = 0;
smname = 'frontrearfun';
modelsb = models;


%parameter optimization

%Levenberg-Marquardt
info.f = smname;
info.tr = 1;
info.tolf = 1e-4;
info.tolx = sqrt(info.tolf);
info.maxit = 300;
info.nu0 = .01;
info.jac = 0;
info.prt = 2;
tic
[x, J, ff, g, iter, conf] = marqdt(info, xv, y, xf, chb, models);
toc
xx = x0;
xx(pvar) = x; %estimated parameters

%
% get residuals and estimator
%
[F, e, hb, Mb, Pevf, A, P] = eval([smname, '(x,y,xf,1,models);']);

% compute residual diagnostics
Ss = e' * e;
Ff = F' * F;
ne = length(e); %residual sum of squares
% disp('concentrated parameter:')
conp = Ss / (ne - nr - nbeta); %estimated sigma square
% disp('square root of concentrated parameter:')
sconp = sqrt(conp);
% compute prediction error variance (finite sample)
% disp('prediction error variance (finite sample)')
Pevf = Pevf * conp;
% disp('standard error (finite sample)')
SPevf = sqrt(diag(Pevf));
%
% lagl=min(36,max([floor(.2*ny) 3*s 10]));
lagl = 36;
ndrs = ne + nbeta;
infr = rescomp(e, lagl, nr, Ss, conp, sconp, Ff, ndrs, nbeta);


%standard errors via second derivatives of log likelihood
xt = x';
H = fdhess('logF', xt, smname, y, xf, 0, models);
SS = inv(H/2) / (ne - nr);
se = sqrt(abs(diag(SS)))';
%t-values
% disp('t-values:')
tt = zeros(size(xx));
tt(pfix) = NaN;
tt(pvar) = x ./ se;
%regression part
if (nbeta > 0)
    seb = sqrt(diag(Mb*conp));
    tb = hb ./ seb; %standard errors and t-values
else
    seb = [];
    tb = [];
end


disp(' ');
disp('******************** Results from estimation ********************');
disp(' ');
result.xvf = xx(pvar);
result.xf = xx(pfix);
result.e = e;
result.Ss = Ss;
result.Ff = Ff;
result.sigma2c = conp;
result.Pevf = Pevf;
result.SPevf = SPevf;
result.tv = tt(pvar)';
result.se = se;
result.F = F;
result.h = hb;
result.M = Mb;
result.A = A;
result.P = P;
result.tvr = tb;
resultser = seb;
%      ferror: 0result.tv=tt';
mprintr(result)
disp(' ')
fprintf(1, '%s %9.4f\n', 'Concentrated parameter:', sconp);
disp('(concentrated parameter is a standard deviation)')
disp('press any key to continue')
pause

disp('More estimation and diagnostic details are in file "frontrear.txt"')
disp('in the subdirectory "results"')
disp('press any key to continue')
pause


%file for output
fid = fopen(fname, 'w');
% fid=1;
%print estimation results
inft = minft(fid, 1, 9, 3, 1);
vnames = char('front','rear');
ny = size(y, 1);
lam = 1;
prtser(fid, vnames(1, :), y1, [], ny, dateim, inft, lam);
prtser(fid, vnames(2, :), y2, [], ny, dateim, inft, lam);

fprintf(fid, 'Estimation results:\n');
stord = models.stord;
conc = models.conc;
arp = 0;
nst = length(stord) - arp;
if nst > 0, xx(1:end-arp) = xx(1:end-arp) * sconp;
end


sse = zeros(size(xx));
sse(pfix) = NaN;
sse(pvar) = se;
xx = [sconp, xx];
sse = [NaN, sse];
tt = [NaN, tt];

disp('Estimated covariance matrix of trend')
S = [xx(1), 0; xx(2), xx(3)];
S * S'
disp('press any key to continue')
pause
disp('Estimated covariance matrix of seasonal')
T = [xx(4), 0; xx(5), xx(6)];
T * T'
disp('press any key to continue')
pause

disp('Estimated covariance matrix of irregular')
M = [xx(7), 0; xx(8), xx(9)];
M * M'
disp('press any key to continue')
pause


z = [xx', sse', tt'];
clear in
in.cnames = char('  Estimate', 'Std. Error', '   T-ratio');
in.rnames = ['Parameter      '];
rnamess = ['Sigma level 11 '; 'Sigma level 21 '; ...
    'Sigma level 22 '; 'Sigma seas  11 '; ...
    'Sigma seas  21 '; 'Sigma seas  22 '; ...
    'Sigma irreg 11 '; 'Sigma irreg 21 '; ...
    'Sigma irreg 22 '; ...
    ];
in.rnames = [in.rnames; rnamess(conc, :)];
if nst > 0, for i = 1:nst, in.rnames = [in.rnames; rnamess(stord(i), :)];
    end, end
in.fmt = char('%12.4f', '%12.4f', '%12.4f');
in.fid = fid;
mprint(z, in);
fprintf(fid, '\nResidual standard error:%11.4f', SPevf);
fprintf(fid, ['\nParameter ', rnamess(conc, :), ' is concentrated out of the likelihood\n']);

%print regression variables
if nbeta > 0
    fprintf(fid, '\n');
    fprintf(fid, 'Regression parameters:\n');
    Mbeta = [hb, seb, tb];
    clear in
    in.cnames = char('  Estimate', 'Std. Error', '   T-ratio');
    rnames = char('Parameter');
    for i = 1:nbeta, rnames = char(rnames, ['reg', num2str(i)]);
    end
    in.rnames = rnames;
    in.fmt = char('%12.5f', '%12.5f', '%8.2f');
    in.fid = fid;
    mprint(Mbeta, in);
    fprintf(fid, '\n');
end
%print residual diagnostics
printres(fid, infr);
%end of residual diagnostics


%plot residual diagnostics
plotres([], [], [], [], [], 1.96, 'residuals', 1, 0, [], 0, [], infr, 1, 1);
pause
close all
