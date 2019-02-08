%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example 6.1 of Tsay (2014), pp. 335-341, with some missing observations
%
% Monthly housing data of the United States from
% January 1963 to July 2012. The two series employed are
% 1. z1t: Logarithm of new homes sold in thousands of units (new residential sales)
% 2. z2t: Logarithm of the total new privately owned housing units started in
%             thousands of units (new residential construction)
% The aim of this exercise is to estimate a multivariate ARMA(p,q)(P,Q)
% with missing observations. To this end, we replace first the missing
% observations with tentative values. In this way, we can identify and
% estimate a model using the Hannan-Rissanen method. Later, in a second
% stage, we can estimate the model with the missing values using the Kalman
% filter. We can also interpolate the missing values using the KF.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

da = load(fullfile('data', 'm-hsoldhst6312.dat'));
zto = log(da(:, 3:4));
tdx = da(:, 1) + da(:, 2) / 12;
[nx, mx] = size(zto);
freq = 12;

%make some observations missing in each series
ztm = zto;
ztm(2, 1) = NaN;
ztm(34, 1) = NaN;
ztm(54:60, 1) = NaN(7, 1);
ztm(350, 1) = NaN;
ztm(401, 1) = NaN;
ztm(503, 1) = NaN;
ztm(2, 2) = NaN;
ztm(37, 2) = NaN;
ztm(154:160, 2) = NaN(7, 1);
ztm(250, 2) = NaN;
ztm(301, 2) = NaN;
ztm(503, 2) = NaN;
ztm(514, 2) = NaN;

%fill in the missing observations with tentative values
[zt1, Xm1, nmiss1, idxn1] = chmarima(ztm(:, 1));
[zt2, Xm2, nmiss2, idxn2] = chmarima(ztm(:, 2));
nmiss = nmiss1 + nmiss2;
zt = [zt1, zt2]; %filled in series
idx = find(~ismember(idxn2, idxn1));
x = [Xm1, Xm2(:, idx)]; %indicator variables for complete missing observations


subplot(2, 1, 1)
plot(tdx, zt(:, 1))
xlabel('time');
ylabel('filled in Hsold');
axis('tight');
subplot(2, 1, 2)
plot(tdx, zt(:, 2))
xlabel('time');
ylabel('filled in Hstart');
axis('tight');
disp('press any key to continue')
pause
close all

%difference all series
yd = diferm(zt, freq); %seasonal differencing
yd = diferm(yd, 1); %regular differencing
xd = diferm(x, freq); %seasonal differencing
xd = diferm(xd, 1); %regular differencing

% %we approximately estimate by regression the missing values. Actually, it
% %should be a restricted regression because there are partial missing
% %values, but the approximation is good enough.
% ct=ones(size(xd,1),1);
% beta=mulols(yd,[xd ct]);
% %and use the corrected differenced series as new filled in series
% yd=yd - [xd ct]*beta;

disp(' ')
disp('estimate VARMA(0,3)(0,1) model using the Hannan-Rissanen method ')
disp('press any key to continue')
pause

%estimate a VARMAX(0,3,0)(0,1,0) model by the Hannan-Rissanen method.
x = [];
hr3 = 0;
finv2 = 1;
mstainv = 1;
nsig = [1, 1];
tsig = [.8, .8];
qr = 3;
qs = 1;
[strv, ferror] = estvarmaxpqrPQR(yd, x, freq, [0, qr, 0], [0, qs, 0], hr3, finv2, mstainv, nsig, tsig);

disp(' ');
disp('***** Estimated VARMA(0,3)(0,1) Model using the HR method  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'th';
strt = 1;
mprintar(strv.thetas3(:, :, 2:qr+1), in, tit, strt);
disp(' ');
tit = 'Th';
strt = 1;
mprintar(strv.thetas3(:, :, freq+qs), in, tit, strt);
disp('press any key to continue')
pause

%estimate using exact ML with the series containing missing observations
%setup model
Phi = eye(2);
Th(:, :, 1) = eye(2);
Th(:, :, qs+1) = strv.thetas3(:, :, freq+qs);
th = strv.thetas3(:, :, 1:qr+1);
phi = eye(2);
Sigma = strv.sigmar3;

[str, ferror] = suvarmapqPQ(phi, th, Phi, Th, Sigma, freq);

%eliminate insignificant parameters
for k = 2:qr + 1
    for i = 1:mx
        for j = 1:mx
            if th(i, j, k) == 0
                str.thn(i, j, k) = 0;
            end
        end
    end
end
for i = 1:mx
    for j = 1:mx
        if Th(i, j, 2) == 0
            str.Thn(i, j, 2) = 0;
        end
    end
end
[str, ferror] = fixvarmapqPQ(str);

Y = [];

disp(' ')
disp('estimation using the exact method with missing observations')
disp('press any key to continue')
pause

ydm = diferm(ztm, freq); %seasonal differencing
ydm = diferm(ydm, 1); %regular differencing
str.nmiss = 1;

%estimate model using the exact method
result = varmapqPQestim(ydm, str, Y);

disp(' ');
disp('******************** Results from estimation ********************');
disp(' ');
mprintr(result)
disp('press any key to continue')
pause

%estimated and fixed parameters
xvf = result.xvf;
xf = result.xf;
%t-values of varma estimated parameters are in result.tv

%create estimated model
[phif, thf, Phif, Thf, Lf, ferror] = pr2varmapqPQ(xvf, xf, str);
Sigmar = result.Sigmar;
%t-values
tvf = result.tv;


[phitvf, thtvf, Phitvf, Thtvf, Ltvf, ferror] = pr2varmapqPQ(tvf, xf, str);

disp(' ');
disp('***** Estimated Model With Missing Observations *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'th';
strt = 1;
mprintar(thf(:, :, 2:qr+1), in, tit, strt);
disp(' ')
tit = 'Th';
strt = 1;
mprintar(Thf(:, :, 2), in, tit, strt);
tit = 'Sigma:';
mprintar(Sigmar, in, tit);
disp(' ')
disp('press any key to continue')
pause

disp(' ');
disp('***** t-values  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'tv-th';
strt = 1;
mprintar(thtvf(:, :, 2:qr+1), in, tit, strt);
disp(' ')
tit = 'tv-Th';
strt = 1;
mprintar(Thtvf(:, :, 2), in, tit, strt);
tit = 'tv-L';
mprintar(Ltvf, in, tit);
disp('press any key to continue')
pause

disp(' ')
disp('******** Computation of interpolated values   ********');
disp(' ');
disp('press any key to continue')
pause

%set up state space model for the original series
dr(:, :, 1) = eye(mx);
dr(:, :, 2) = -eye(mx);
phi = pmatmul(phif, dr);
Phi = pmatmul(Phif, dr);
%pass model to state space form
[~, ~, Z, T, H, G, ~] = varmapqPQ2ssm(phi, thf, Phi, Thf, Lf, str);
%
% Computation with function smoothgen.m
%
% Function smoothgen smooths a general vector:
% Y_t = U_t*beta + C_t*alpha_t + D_t*epsilon_t
% In this case, it is desired to smooth:
% Y_t = Z*alpha_t + G*epsilon_t
% Hence, U_t = 0, C_t=Z and D_t=G
U = [];
mucd = mx;
C = Z;
D = G;
ndelta = (1 + freq) * mx;
[ins, ii, ferror] = incossm(T, H, ndelta);
[ztmint, sztmint] = smoothgen(ztm, [], Z, G, [], T, H, ins, ii, mucd, U, C, D);

%obtain interpolated and original values
Interp1 = ztmint(idxn1, 1);
Orig1 = zto(idxn1, 1);
Interp2 = ztmint(idxn2, 2);
Orig2 = zto(idxn2, 2);
sInterp1 = Interp1;
for i = 1:nmiss1
    sInterp1(i) = sztmint(idxn1(i)*mx-1, 1);
end
sInterp2 = Interp2;
for i = 1:nmiss2
    sInterp2(i) = sztmint(idxn2(i)*mx, 2);
end
sInterp1 = sqrt(sInterp1*result.sigma2c);
sInterp2 = sqrt(sInterp2*result.sigma2c);
disp(' ');
disp('***** Interpolated values  *****');
disp(' ');
clear in
in.cnames = char('  Estimate', 'Std. Error', '  Original value');
rnamesrgi = ['interp1. ', num2str(idxn1(1))];
for i = 2:nmiss1
    rnamesrgi = char(rnamesrgi, ['interp1. ', num2str(idxn1(i))]);
end
rnames = char('Interpolated value (series 1)', rnamesrgi);
in.rnames = rnames;
in.fmt = char('%12.5f');
mprint([Interp1, sInterp1, Orig1], in);
disp(' ')
clear in
in.cnames = char('  Estimate', 'Std. Error', '  Original value');
rnamesrgi = ['interp2. ', num2str(idxn2(1))];
for i = 2:nmiss2
    rnamesrgi = char(rnamesrgi, ['interp2. ', num2str(idxn2(i))]);
end
rnames = char('Interpolated value (series 2)', rnamesrgi);
in.rnames = rnames;
in.fmt = char('%12.5f');
mprint([Interp2, sInterp2, Orig2], in);
disp('press any key to continue')
pause

disp(' ')
disp('******** Computation of smooth Residuals   ********');
% In this case, it is desired to smooth:
% Y_t =  G*epsilon_t
nalpha = size(T, 1);
C = zeros(mx, nalpha);
D = G;
[sres, ssres] = smoothgen(ztm, [], Z, G, [], T, H, ins, ii, mucd, U, C, D);
%plot smooth residuals
plot(sres(:, 1)), legend('Hsold'), pause
plot(sres(:, 2)), legend('Hstart'), pause
close all
%compute autocovariance and autocorrelation matrices of smooth residuals
lag = 24;
ic = 1;
nr = 0; %length(result.xvf);
disp(' ')
disp('******** Smooth Residuals:     ********');
str = mautcov(sres, lag, ic, nr);
disp('Correlation matrix at lag 0:')
disp(str.r0)
disp('Q statistics:')
disp(str.qstat)

disp('p-values of Q statistics:')
disp(str.pval)
[m, n] = size(str.pval);
t = 1:m;
plot(t, str.pval, t, 0.05*ones(1, m))
legend('p-values of Q statistics:')
pause
close all
