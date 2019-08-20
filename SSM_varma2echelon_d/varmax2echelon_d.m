%
%script file to illustrate the computation of the echelon form of a VARMAX
%model or a state space model
%
%1) VARMAX model
% Kronecker indices [1 1 0]
%
clear
phi0 = [1., 0., 0.; 0., 1., 0.; -.1, .2, 1.];
phi1 = [.7, 0., 0.; .4, .9, 0.; 0., 0., 0.];
th1 = [.5, .4, -.1; 0., .5, .3; 0., 0., 0.];
g0 = [-2., 1.; .5, 3.; -.5, 1.];
g1 = [3., -1.; 2., 1.; 0., 0.];
phi = cat(3, phi0, phi1);
theta = cat(3, phi0, th1);
gamma = cat(3, g0, g1);
disp(' ');
disp('***** Original Model  *****');
disp(' ');
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'phi';
strt = 1;
mprintar(phi, in, tit, strt);
disp(' ')
tit = 'theta';
mprintar(theta, in, tit, strt);
disp(' ')
tit = 'gamma';
mprintar(gamma, in, tit, strt);
disp('press any key to continue')
pause
Phi(:, :, 1) = eye(3);
Phi(:, :, 2) = phi0 \ phi1;
Th(:, :, 1) = eye(3);
Th(:, :, 2) = phi0 \ th1;
Gamma(:, :, 1) = phi0 \ g0;
Gamma(:, :, 2) = phi0 \ g1;
disp(' ');
disp('***** Model premultiplied by phi0^{-1}  *****');
disp(' ');
tit = 'Phi';
strt = 1;
mprintar(Phi, in, tit, strt);
disp(' ')
tit = 'Th';
mprintar(Th, in, tit, strt);
disp(' ')
tit = 'Gamma';
mprintar(Gamma, in, tit, strt);
disp('press any key to continue')
pause
Theta(:, :, 1) = [Gamma(:, :, 1), Th(:, :, 1)];
Theta(:, :, 2) = [Gamma(:, :, 2), Th(:, :, 2)];
[phie, thetae, kro, ierror] = pecheform(Phi, Theta);
disp(' ');
disp('***** Model in Echelon Form *****');
disp(' ');
tit = 'kro';
mprintar(kro, in, tit);
disp(' ')
tit = 'phie';
strt = 1;
mprintar(phie, in, tit, strt);
disp(' ')
tit = 'thetae';
mprintar(thetae, in, tit, strt);
disp('press any key to continue')
pause

%
%2) state space model
%
%define univariate structural model: trend, slope, trigonometric
%seasonality, cycle, irregular and autoregressive component
comp.level = [1, 1., NaN];
comp.slope = [1, 0.2, NaN];
comp.seas = [2, .5, NaN];
comp.irreg = [1, .3, NaN];
freq = 4;
comp.freq = freq;
bg_year = 1981;
bg_per = 1;
datei = cal(bg_year, bg_per, freq);
comp.datei = datei;
y = [];
npr = 0;
Y = [];

%create structure and put model into state space form
[str, ferror] = suusm(comp, y, Y, npr);
disp(' ');
disp('***** Structural Model *****');
disp(' ');
tit = 'T';
mprintar(str.T, in, tit);
disp(' ')
tit = 'H';
mprintar(str.H, in, tit);
disp(' ')
tit = 'Z';
mprintar(str.Z, in, tit);
disp(' ')
tit = 'G';
mprintar(str.G, in, tit);
disp('press any key to continue')
pause
%reduced form of the estimated structural model
T = str.T;
G = str.G;
Z = str.Z;
H = str.H;
nalpha = size(T, 1);
[ng, mg] = size(G);
phir(:, :, 1) = eye(nalpha);
phir(:, :, 2) = -T;
thetar(:, :, 2) = Z;
np = nalpha + 1;
[phie, thetae, kro, ierror] = pright2leftcmfd(phir, thetar, np);
th = pmatmul(phie, G) + pmatmul(thetae, H);
cov = pmmulbf(th, th);
[nc, mc, pc] = size(cov);
ninit = floor(pc/2) + 1;
Lp = cov(:, :, ninit:end);
[Omega, Theta, ierror, iter, normdif] = pmspectfac(Lp, 30);
%phie contains the AR part, Theta contains the MA part except Theta(0)=1
thetaec = phie;
thetaec(:, :, 2:end) = Theta;
disp(' ');
disp('***** Model in Echelon Form *****');
disp(' ');
tit = 'kro';
mprintar(kro, in, tit);
tit = 'phie';
strt = 1;
mprintar(phie, in, tit, strt);
disp(' ')
tit = 'thetaec';
mprintar(thetaec, in, tit, strt);
disp(' ')
tit = 'Omega';
mprintar(Omega, in, tit);
disp('press any key to continue')
pause
%
%alternative way to compute the innovations form using the DARE
%
%pass structural model to innovations form
%
[P, K, Sigma, U, iU] = ss2if(T, H, Z, G, eye(mg), zeros(mg), eye(mg));
%obtain moving average polynomial. AR polynomial is the same
thdare = phie + pmatmul(thetae, K);
disp(' ');
disp('***** Model in Echelon Form Using the DARE *****');
disp(' ');
tit = 'phie';
strt = 1;
mprintar(phie, in, tit, strt);
disp(' ')
tit = 'thdare';
mprintar(thdare, in, tit, strt);
%innovations variance computed with the DARE
disp(' ')
tit = 'Sigma';
mprintar(Sigma, in, tit);
