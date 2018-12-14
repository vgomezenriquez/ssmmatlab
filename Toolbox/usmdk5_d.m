%Example of a spline smoothing model.
%Series is car drivers killed or seriously injured in Great Britain from
%January 1969 to December 1984 (Durbin and Koopman, 2012)
%
%model is
%
%         y_i  = \mu_i + \epsilon_i
%
% [\mu_{i+1}]  = [1  \delta_i][\nu_i]   + [\xi_i ]
% [\nu_{i+1}]    [0      1   ][\nu_i]     [\eta_i]
%
% Var(\epsilon_i)=\sigma^2_\epsilon
%
%  Q_i = \sigma^2_\eta*\delta_i[\delta^2_i/3    \delta_i/2]
%                              [\delta_i/2           1    ]
%
%  \lambda= \sigma_\eta/\sigma_\epsilon
%
% concentrated parameter: \sigma_\epsilon
%
% I believe what is reported in (Durbin and Koopman, 2012) is \lambda
%(0.0275) instead of the smoothing parameter. Indeed, in Harvey and Koopman
% (2000), pp. 98 and 100, the parameter \gamma, which is the equivalent to
% \lambda here, is the quotient of standard deviations.
%

clear

yy = load(fullfile('data', 'motorcycle.dat'))';
ss = yy(1, :);
dt = yy(2, [2:end, 1])';
y = yy(3, :)';
yor = y;
tname = 'motorcycle';
fname = fullfile('results', 'motorcycle.txt');
lam = 1;
npr = 12; %number of forecasts

lag = 36;
cw = 1.96;
freq = 1;
dr = 0;
ds = 0;
c0 = sacspacdif(y, tname, dr, ds, freq, lag, cw);
pause
close all

%set up matrices for state space model
ny = size(y, 1);
ns = 2;
ne = 3;
T = zeros(ns*ny, ns);
H = zeros(ns*ny, ne);
for i = 1:ny
    %  i
    ip = (i - 1) * ns + 1:i * ns;
    T(ip, :) = [1, dt(i); 0, 1];
    if dt(i) ~= 0
        M = dt(i) * [dt(i)^2 / 3, dt(i) / 2; dt(i) / 2, 1];
        [U, D, V] = svd(M);
        H(ip, 1:ns) = U * sqrt(D) * U';
    end
end
Z = [1, 0];

x0 = .1;
xv = x0;

G = [0, 0, 1.];
Y = [];
W = [];
ib = [0, 0, 0, 2];
insb = eye(2);
models.X = Y;
models.W = W;
models.T = T;
models.H = H;
models.Z = Z;
models.G = G;
models.i = ib;
models.ins = insb;
nbeta = 0;
nr = 1;

clear T H Z G ib insb;


s = freq;
chb = 0;
smname = 'motorcyclefun';


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
[x, J, ff, g, iter, conf] = marqdt(info, xv, y, chb, models);
toc
xx = x;


%
% get residuals and estimator
%
[F, e, hb, Mb, Pevf, A, P] = eval([smname, '(x,y,1,models);']);

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
H = fdhess('logF', xt, smname, y, 0, models);
SS = inv(H/2) / (ne - nr);
se = sqrt(abs(diag(SS)))';
%t-values
% disp('t-values:')
tt = x ./ se;
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
result.xvf = xx;
result.xf = [];
result.e = e;
result.Ss = Ss;
result.Ff = Ff;
result.sigma2c = conp;
result.Pevf = Pevf;
result.SPevf = SPevf;
result.tv = tt;
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
disp('press any key to continue')
pause

disp('More diagnostic details are in file "motorcycle.txt"')
disp('in the subdirectory "results"')
disp('press any key to continue')
pause


%plot residual diagnostics
plotres([], [], [], [], [], 1.96, 'residuals', 1, 0, [], 0, [], infr, 1, 1);
pause
close all

%file for output
fid = fopen(fname, 'w');
%print residual diagnostics
printres(fid, infr);
%close external file
if fid ~= 1
    fclose(fid);
end
pause
close all


sigma2eps = conp;
sigma2eta = xx^2 * conp;
disp('standard deviation ratio:')
lambda = abs(xx)
disp('smoothing paramter:')
lambda^2


%smoothing
Z = models.Z;
T = models.T;
G = models.G;
H = models.H;
W = models.W;
X = models.X;
ins = models.ins;
ii = models.i;

H = H * xx;

[Xt, Pt, g, M] = scakfs(y, X, Z, G, W, T, H, ins, ii);

cw = 1.69;
ny = size(y, 1);
trend = Xt(:, 1);
strend = zeros(size(trend));
for i = 1:ny
    strend(i) = sqrt(Pt((i - 1)*ns+1, 1)) * sconp;
end

figure
plot(ss, y, ss, zeros(length(ss)), ss, trend, ss, trend+cw*strend, ':', ...
    ss, trend-cw*strend, ':')
legend('Original series with trend and 95% confidence interval')

disp('strike a key to continue')
pause
close all
