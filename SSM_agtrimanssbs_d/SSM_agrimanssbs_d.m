%
% Practical example for the paper ''A new State Space Methodology to
% Disaggregate Multivariate Time Series'', Journal of Time Series Analysis,
% 30, 97-124, by Gomez and Aparicio, (2009).
% The basic solution is specified by setting theta=0
% The specific solution is specified by setting theta non equal to zero.
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
nbeta = 0;

ipi_19951_200712 = load(fullfile('data', 'ipi_19951_200712.dat'));
pib_1995I_2007IV = load(fullfile('data', 'pib_1995I_2007IV.dat'));
ipi = ipi_19951_200712(:, 2);
pib = pib_1995I_2007IV(:, 2);
p = 4;
lam = 1;


%plot both series
bg_yearm = 1995;
bg_perm = 1;
freqm = 12;
dateim = cal(bg_yearm, bg_perm, freqm);
tsplot(ipi, dateim, 'IPI');
disp('strike a key to continue')
pause
bg_yeart = 1995;
bg_pert = 1;
freqt = 4;
dateit = cal(bg_yeart, bg_pert, freqt);
tsplot(pib, dateit, 'GDP');
disp('strike a key to continue')
pause

%compute IPI quarterly figures
[nq, junk] = size(pib);
ipiq = zeros(size(pib));
for i = 1:nq
    ipiq(i) = mean(ipi((i - 1)*3+1:i*3));
end
tsplot(ipiq, dateit, 'IPIQ');
disp('strike a key to continue')
pause

yq = [pib, ipiq];
vnames =char('GDP','Quarterly aggregated IPI');
tsplot(yq, dateit, vnames);
disp('strike a key to continue');
pause

%compute GDP annual figures
ny = floor(nq/4);
piby = zeros(ny, 1);
for i = 1:ny
    piby(i) = mean(pib((i - 1)*4+1:i*4));
end
bg_yeary = 1995;
bg_pery = 1;
freqy = 1;
dateiy = cal(bg_yeary, bg_pery, freqy);
tsplot(piby, dateiy, 'GDPY');
disp('strike a key to continue')
pause

nq = 4 * ny;
ipiq = ipiq(1:nq, :);

%compute bivariate series
yq = [];
for i = 1:ny
    yq = [yq; [ipiq(4*(i - 1)+1:4*i)', piby(i)]];
end

tname = 'agtrimanss1';
fname = fullfile('results', 'agtrimanss1.txt');

%first, form univariate structural model and, from this, construct
%bivariate structural model (this is easier than constructing the bivariate
%structural model from scratch)

%univariate structural model defined according to ssm_matlab
% codes for the components of the univariate structural model:
% level = -1  constant
%          1  stochastic structural (Butterworth sine)
%          2  Butterworth tan
% slope = -1  constant
%          0  zero
%          1  stochastic structural (Butterworth sine)
% seas  = -1  fixed dummy seasonality
%          1  stochastic dummy seasonality
%          2  trigonometric seasonality
%          4  Butterworth tan seasonality
% irreg =  1  stochastic
% arp   =  k   autoregressive of order k
% (see ssm_matlab manual for further information)
%
Y = [];
comp.level = [1, .1, NaN];
comp.slope = [1, .1, NaN];
comp.seas = [2, .1, NaN];
comp.irreg = [1, .1, NaN];
comp.conout = 'level';
comp.freq = freqt;
comp.datei = dateit;
npr = 0;
[models, ferror] = suusm(comp, ipiq, Y, npr);
%end of univariate structural model defined according to ssm_matlab
nreg = 0;
nr = length(models.pvar);

% Bivariate structural model
Tb = kron(models.T, eye(2));
Hb = kron(models.H, eye(2));
Zb = kron(models.Z, eye(2));
Gb = kron(models.G, eye(2));
insb = kron(models.ins, eye(2));
[nt, mt] = size(models.T);
ib = [0, 0, 0, 2 * nt];

% order in the standard deviations (conc):
% 1 2 3    - level
% 4 5 6    - slope
% 7 8 9    - seasonal
% 10 11 12 - irregular
x0 = [.1, .1, ... %sigma level s21 s22 (s11 one parameter is concentrated out)
    .1, .1, .1, ... %sigma slope s11 s21 s22
    .1, .1, 0, ... %sigma seasonal s11 s21 s22 (s22 is not identified)
    .1, .1, .1, ... %sigma irreg. s11 s21 s22
    ];

pfix = [8]; %fixed parameters
pvar = [1:7, 9:11]; %free parameters

%parameter for basic solution
% theta=0.0;

%parameter for specific solution
% theta=0.5271;  %parameter for proportional seasonality (0.2002/1.8767+0.0361/0.0381)/2
theta = 0.5903; %parameter for proportional seasonality (0.1280/1.4271+0.0024/0.0022)/2

% theta=.75232;    %parameter for proportional seasonality, Fernandez

xv = x0(pvar);
xf = x0(pfix);

stordb = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]; %standard deviations to be estimated
concb = 1; %standard deviation concentrated out
[nhb, mhb] = size(Hb);
xp = zeros(1, 12);
xp(stordb) = x0;
xp(concb) = 1; %concentrated parameter
Hb(1, 1) = xp(1);
Hb(2, 1) = xp(2);
Hb(2, 2) = xp(3);
Hb(3, 3) = xp(4);
Hb(4, 3) = xp(5);
Hb(4, 4) = xp(6);
j = 5;
for i = 6:2:nhb
    Hb(i, j) = xp(8);
    Hb(i-1, j) = xp(7);
    Hb(i, j+1) = xp(9);
    j = j + 2;
end
% for i=6:2:nhb
%  Hb(i,5)=xp(8); Hb(i-1,5)=xp(7); Hb(i,6)=xp(9);
% end
Gb(1, mhb-1) = xp(10);
Gb(2, mhb-1) = xp(11);
Gb(2, mhb) = xp(12);

models.T = Tb;
models.H = Hb;
models.Z = Zb;
models.G = Gb;
models.i = ib;
models.ins = insb;
models.stord = stordb;
models.conc = concb;
% end of bivariate structural model

clear Tb Hb Zb Gb ib insb stordb concb;

%aggregation pattern
D = [1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 1, 0, 0, 0, 0, 0; 0, 0, 0, 0, 1, 0, 0, 0; ...
    0, 0, 0, 0, 0, 0, 1, 0; 0, 1 / 4, 0, 1 / 4, 0, 1 / 4, 0, 1 / 4];
models.D = D;

s = freqt;
stcs = 1;
chb = 0;
Q = [];
models.nalphao = [];
[modelsn, H_p, J_p, Q] = modstr_agtrimanbs(xv, pfix, pvar, xf, models, Q, stcs, p, ny);
models.nalphao = modelsn.nalphao;
Z = modelsn.Z;
T = modelsn.T;
G = modelsn.G;
H = modelsn.H;
W = modelsn.W;
X = modelsn.X;
ins = modelsn.ins;
i = modelsn.i;
[nalphao, junk] = size(T);
smname = 'agtrimanbsfun';
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
[x, J, ff, g, iter, conf] = marqdt(info, xv, yq, pfix, pvar, xf, 0, modelsb, Q, p, ny);
toc
xx = x0;
xx(pvar) = x; %estimated parameters


%
% get residuals and estimator
%
[F, e, hb, Mb, Pevf, A, P] = eval([smname, '(x,yq,pfix,pvar,xf,1,modelsb,Q,p,ny);']);

% compute residual diagnostics
Ss = e' * e;
Ff = F' * F;
ne = length(e); %residual sum of squares
disp('concentrated parameter:')
conp = Ss / (ne - nr); %estimated sigma square
disp('square root of concentrated parameter:')
sconp = sqrt(conp);
% compute prediction error variance (finite sample)
disp('prediction error variance (finite sample)')
Pevf = Pevf * conp;
disp('standard error (finite sample)')
SPevf = sqrt(diag(Pevf));
%
% lagl=min(36,max([floor(.2*ny) 3*s 10]));
lagl = 36;
ndrs = ne + nbeta;
infr = rescomp(e, lagl, nr, Ss, conp, sconp, Ff, ndrs, nbeta);

%standard errors via second derivatives of log likelihood
xt = x';
H = fdhess('logF', xt, smname, yq, pfix, pvar, xf, 0, modelsb, Q, p, ny);
SS = inv(H/2) / (ne - nr);
se = sqrt(abs(diag(SS)))';
%t-values
disp('t-values:')
tt = zeros(size(xx));
tt(pfix) = NaN;
% tt(pvar)=x./se

%file for output
fid = fopen(fname, 'w');
% fid=1;
%print estimation results
inft = minft(fid, 1, 9, 3, 1);
vnames = char('IPI','GDP');
prtser(fid, vnames(1, :), ipiq, [], nq, dateit, inft, lam);
prtser(fid, vnames(2, :), piby, [], ny, dateiy, inft, lam);

fprintf(fid, 'Estimation results:\n');
stord = models.stord;
conc = models.conc;
arp = models.arp;
nst = length(stord) - arp;
if nst > 0, xx(1:end-arp) = xx(1:end-arp) * sconp;
end
sse = zeros(size(xx));
sse(pfix) = NaN;
sse(pvar) = se;
xx = [sconp, xx];
sse = [NaN, sse];
tt = [NaN, tt];
z = [xx', sse', tt'];
clear in
in.cnames = strvcat('  Estimate', 'Std. Error', '   T-ratio');
in.rnames = ['Parameter      '];
rnamess = ['Sigma level 11 '; 'Sigma level 21 '; ...
    'Sigma level 22 '; 'Sigma slope 11 '; ...
    'Sigma slope 21 '; 'Sigma slope 22 '; ...
    'Sigma seas. 11 '; 'Sigma seas. 21 '; ...
    'Sigma seas. 22 '; 'Sigma irreg 11 '; ...
    'Sigma irreg 21 '; 'Sigma irreg 22 '; ...
    ];
in.rnames = [in.rnames; rnamess(conc, :)];
if nst > 0, for i = 1:nst, in.rnames = [in.rnames; rnamess(stord(i), :)];
    end, end
in.fmt = strvcat('%12.4f', '%12.4f', '%12.4f');
in.fid = fid;
mprint(z, in);
fprintf(fid, '\nResidual standard error:%11.4f', SPevf);
fprintf(fid, ['\nParameter ', rnamess(conc, :), ' is concentrated out of the likelihood\n']);

%print residual diagnostics
printres(fid, infr);
%end of residual diagnostics


%plot residual diagnostics
plotres([], [], [], [], [], 1.96, 'residuals', 1, 0, [], 0, [], infr, 1, 1);


%smoothing starts
% [Xt,Pt,g,M]=eval(['agtrimanbssmt' '(x,yq,s,pfix,pvar,xf,1,modelsb,Q,p,ny);']); %estados
% iti=11; on=2;
% format long;
% Xt(on,:)
% conp*Pt((on-1)*nalphao+1:on*nalphao,:)
% format short;


%interpolation. Update system matrices
stcs = 0;
[modelsn, H_p, J_p, Q] = modstr_agtrimanbs(x, pfix, pvar, xf, modelsb, Q, stcs, p, ny);
Ca = H_p * Q';
D = J_p;
if ~isempty(Y)
    U = ypr;
else
    U = [];
end
C = Ca(:, end-nalphao+1:end);
[mucd, junk] = size(C);

%specific solution

beta = [0, 0, 0, 0, theta, 0, 0; 0, 0, 0, 0, 0, theta, 0; 0, 0, 0, 0, 0, 0, theta];

Qs = Q(end-nalphao+1:end, [1:5, 7, 9]) + Q(end-nalphao+1:end, [6, 8, 10]) * beta;
Cn = Q(1:end-nalphao, [1:5, 7, 9]) + Q(1:end-nalphao, [6, 8, 10]) * beta;
Cn = Cn / Qs;
C = C + Ca(:, 1:end-nalphao) * Cn;
%end of specific solution


% [junk,ndd]=size(D);
% mucd=nalphao; C=zeros(nalphao); D=zeros(nalphao,ndd); U=ones(nalphao,1); %state vector
% mucd=ndd; C=zeros(ndd,nalphao); D=eye(ndd); U=zeros(ndd,1);     %innovations vector
% U
% C
% D
% C=sparse(C); D=sparse(D); U=sparse(U);

[Xt, Pt, g, M] = eval(['agtrimanbsint', '(x,yq,pfix,pvar,xf,modelsb,Q,mucd,U,C,D,p,ny);']);
% format long;
% Xt(on,:)
% conp*Pt((on-1)*mucd+1:on*mucd,:)
% format short;

ipiint = [];
pibint = [];
ripiint = [];
rpibint = [];
for on = 1:ny
    yint = Xt(on, :);
    ryint = sconp * sqrt((diag(Pt((on - 1)*mucd+1:on*mucd, :)))); %ojo, valores negativos
    %  ipiint=[ipiint; yint(1); yint(3); yint(5)];
    pibint = [pibint; yint(2); yint(4); yint(6); yint(8)];
    %  ripiint=[ripiint; ryint(1); ryint(3); ryint(5)];
    rpibint = [rpibint; ryint(2); ryint(4); ryint(6); ryint(8)];
end

fprintf(fid, '\n');
inft = minft(fid, 1, 9, 3, 1);
% dateib=dateim; inftb=inft; fnameb=' Interpolated IPI'; nint=length(ipiint);
% prtser(fid,fnameb,ipiint,[],nint,dateib,inftb,lam);
dateib = dateit;
inftb = inft;
fnameb = ' Interpolated GDP';
nint = length(pibint);
prtser(fid, fnameb, pibint, [], nint, dateib, inftb, lam);
fprintf(fid, '\n');
% fnameb=' Root Mean Squared Errors: IPI'; nint=length(ripiint);
% prtser(fid,fnameb,ripiint,[],nint,dateib,inftb,lam);
fnameb = ' Root Mean Squared Errors: GDP';
nint = length(rpibint);
prtser(fid, fnameb, rpibint, [], nint, dateib, inftb, lam);
% error bands
% lbripiint=ipiint-1.96*ripiint;
% ubripiint=ipiint+1.96*ripiint;
lbrpibint = pibint - 1.96 * rpibint;
ubrpibint = pibint + 1.96 * rpibint;

%close external file
if fid ~= 1
    fclose(fid);
end

%%To obtain the current ColorOrder, which may be set during startup, get
%%the property value:
% % M=get(gca,'ColorOrder');
%%after this, to restore the colors to the original values, do the
%%following:
% %set(0,'DefaultAxesColorOrder',M)
%
%%the following sets the default ColorOrder to use only the color black and
%%sets the LineStyleOrder to use solid, dash-dot, dash-dash, and dotted line styles.
% set(0,'DefaultAxesColorOrder',[0 0 0],...
%       'DefaultAxesLineStyleOrder','-|--|-.|:')
% figure
% vnames=['             IPI'
%         'Interpolated IPI'];
% tsplot([ipi  ipiint],dateim,vnames);disp('strike a key to continue');
% pause

% ipiqint=zeros(size(ipiq));      % interpolated quarterly IPI. It should be equal
%                                 % to  IPIQ
% for i=1:nq
%  ipiqint(i)=mean(ipiint((i-1)*3+1:i*3));
% end
figure
vnames = char('GDP','Interpolated GDP');
tsplot([pib(1:nq), pibint], dateit, vnames);
disp('strike a key to continue');
pause


% figure
% vnames=['Interpolated IPI',
%         'Interpolated GDP'];
% tsplot([ipiint pibint],dateim,vnames);disp('strike a key to continue');
% pause


figure
vnames = char('Interpolated GDP','Lower band','Upper band');
tsplot([pibint, lbrpibint, ubrpibint], dateit, vnames);
disp('strike a key to continue');

% mspesm = 100*sqrt(sum(((pibint-pib(1:nq))./pib(1:nq)).^2)./nq)
% apesm = 100*sum(abs((pibint-pib(1:nq))./pib(1:nq)))./nq

% msesm = sqrt(sum(((pibint-pib(1:nq))).^2)./nq)
% aesm = sum(abs((pibint-pib(1:nq))))./nq
pause
close all
