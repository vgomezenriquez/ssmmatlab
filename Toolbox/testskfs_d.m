%
%Script file to test the Kalman filter and smoothing functions of SSMMATLAB
%
clear
clc
% 
%part I) 
%set up VARMA model in echelon form. There may be unit roots. If the VARMA
%model is not in echelon form, we can still use the following setup. For
%example, if the number of variables is s = 3, put kro = (k, k, k), where k
%is the degree of the matrix polynomials, in this case phi = I + phi_1 * z  
%+ ... + phi_k * z^k and theta = I + th_1 * z + ... + th_k * z^k. Note that
%phi and theta must have the the same degrees, but this is no restriction
%because they can always be completed with zeros. Note also that the fast
%square root filter can also be used with stationary models with no missing
%values. 
s = 2;
ndelta = 0;                       %number of unit roots
phis(:,:,1) = eye(s); 
phis(:,:,2) = [-.8 0.; 0. -.5];
thetas(:,:,1) = eye(s);
thetas(:,:,2) = zeros(s);
kro = [1 1];
sigma = .3;  
Sigma = eye(s) * (sigma^2);       %innovations variance
LSigma = eye(s) * sigma;          %square root of innovations variance
%create echelon structure str
m = 0;                            %no inputs, no VARMAX but VARMA
gammas = [];   
[str, ~] = matechelon(kro, s, m);
str.phis = phis;
str.thetas = thetas;
str.gammas = gammas;
str.sigmar2 = Sigma; 
str = armaxe2sse(str);

%Part II)
%set up model in state space form. If the model is VARMA in echelon form,
%we can use the previous structure, str. If not, set up the system matrices
%directly instead of the following. 
T = str.Fs;
Z = str.Hs;
H = str.Ks * LSigma;
G = LSigma;
W = [];


%data and regression matrix
X=[1. 0. 0.; 0. 1. 0.; 0. 0. 1.; 1. 0. 0.; 0. 1. 0.; 0. 0. 1.; 1. 0. 0.; 0. 1. 0.];
y=[-1. 1.; 3. -1.; 1. 3.; .5 -1.];
% y=[-1. NaN; 3. -1.; 1. 3.; .5 -1.];
%number of missing observations
missingt = sum(sum(isnan(y)));

%number of observations
n=size(y,1); 

%set up initial conditions
[ins, ii, ferror] = incossm(T, H, ndelta);
% ii(1)=0;
% ins = ins(:,end);
% G=zeros(2);

%estimate the regression coefficients
chb=1;
%filtering t+1|t
[e2, f2, hd2, Md2, A2, P2, qyy2, R2, olsres2] = scakfle2(y, X, Z, G, W, T, H, ins, ii, chb);

[e1, f1, hd1, Md1, A1, P1, n1] = tskfsribf(y, X, Z, G, W, T, H, ins, ii, chb); 

[e3, f3, hd3, Md3, A3, LP, qyy3, R3, olsres3] = scakflesqrt(y, X, Z, G, W, T, H, ins, ii, chb);

%fast square root filtering (stationary models without missing values only)
%this code is a little long at present. But it can be changed easily. 
if (ndelta == 0) && (missingt == 0)
    tol = eps;
    maxupdt = [];  
    [e, E, rSigmat] = sqrt_ckms(y, X, str, maxupdt, tol);
    recr = zeros(n, s);                    %matrix for OLS residuals
    nbeta = size(X,2);
    fsq = 1;                               %determinantal factor
    fc = 0;
    es = zeros(n*s, 1);                    %vector for ml residuals
    SQT = [];
    for jj = 1:n
        ind = (jj - 1) * s + 1:jj * s;
        V = rSigmat(ind, :);
        es(ind) = V \ e(jj, :)';
        if nbeta > 0
            SQT = [SQT; V \ E((jj - 1)*nbeta+1:jj*nbeta, :)'];
        end
        fsq = fsq * abs(prod(diag(V)));
        [fsq, fc] = updatef(fsq, fc);
    end
    nbsqs = n * s;
    fsq = (fsq^(1 / (nbsqs))) * (2^(fc / (nbsqs)));
    R = [];
    SQT = [SQT, es];
    if nbeta > 0
        [Q, R] = qr(SQT(:, 1:nbeta));
        qy = Q' * SQT(:, nbeta+1);
        esq = qy(nbeta+1:end);
        hb = R(1:nbeta, :) \ qy(1:nbeta);
        Mb = pinv(R(1:nbeta, :));
        Mb = Mb * Mb';
    else
        esq = SQT;
        hb = [];
        Mb = [];
    end
    for jj = 1:n
        ind = (jj - 1) * nbeta + 1:jj * nbeta;  
        recr(jj, :) = e(jj, :) - (E(ind, :)' * hb)';
    end
else
    recr = [];
    hb = [];
    esq = [];
    fsq = [];
end
%end of fast square root filtering 

%smoothing
[KKP1, PT1, hd4, Md4] = scakfs(y, X, Z, G, W, T, H, ins, ii);

[KKP2, PT2, hd5, Md5] = scakfssqrt(y, X, Z, G, W, T, H, ins, ii);

mucd = 2;
% U = [1 0 0; 0 1 0];
% C = zeros(2);
% D = zeros(2);
% the following combination is to interpolate y_t. By changing U, C and D,
% the user can smooth other vectors. 
U = X;              
C = Z;
D = G;
[KKP3, PT3, hd6, Md6] = smoothgen(y, X, Z, G, W, T, H, ins, ii, mucd, U, C, D);

%filtering t|t
[KKP4, PT4, hd7, Md7, initf1, recrs1, recr1, srecr1] = scakff(y, X, Z, G, W, T, H, ins, ii); 

[KKP5, PT5, hd8, Md8, initf2, recrs2, recr2, srecr2] = scakffsqrt(y, X, Z, G, W, T, H, ins, ii); 

[KKP6, PT6, recrs3, recr3, srecr3, t1, A4, P4, KG1] = scakfff(y, X, Z, G, W, T, H, ins, ii, hd1);

[KKP7, PT7, recrs4, recr4, srecr4, t13, A5, P5, KG2] = scakfffsqrt(y, X, Z, G, W, T, H, ins, ii, hd1);

%checking fixing regression parameters in filtering t|t
yy = vec(y') - X*hd1;
yy=reshape(yy,[],n)';
XX=[];
[KKP8, PT8, hd9, Md9, initf3, recrs5, recr5, srecr5] = scakff(yy, XX, Z, G, W, T, H, ins, ii);

[KKP9, PT9, hd10, Md19, initf4, recrs6, recr6, srecr6] = scakffsqrt(yy, XX, Z, G, W, T, H, ins, ii);

%results
%sum of squares
esq'*esq
e1'*e1
e2'*e2
e3'*e3

%determinantal terms
fsq
f1
f2
f3

%regression parameters estimates
hb
hd1
hd2
hd3
hd4
hd5
hd6
hd7
hd8
hd9
hd10

%recursive residuals
recr1
recr2
recr
recr3
recr4
recr5
recr6

%filtered or smoothed state vector
KKP1
KKP2
KKP3
KKP4
KKP5
KKP6
KKP7
KKP8
KKP9



