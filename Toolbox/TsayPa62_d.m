%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example in paragraph 6.2 in Tsay (2014), pp. 341-345
%
% 4-dimensional monthly time series zt = (z1t, . . . , z4t) of U.S. manufacturers
% data on durable goods, where
%
% 1. z1t: New orders (NO)
% 2. z2t: Total inventory (TI)
% 3. z3t: Unfilled orders (UO)
% 4. z4t: Values in shipments (VS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

da = load(fullfile('data', 'm-amdur.dat'));
dur = da(:, 3:6) / 1000;
tdx = da(:, 1) + da(:, 2) / 12;


subplot(2, 2, 1)
plot(tdx, dur(:, 1))
xlabel('time');
ylabel('NO');
axis('tight');
subplot(2, 2, 2)
plot(tdx, dur(:, 3))
xlabel('time');
ylabel('UO');
axis('tight');
subplot(2, 2, 3)
plot(tdx, dur(:, 2))
xlabel('time');
ylabel('TI');
axis('tight');
subplot(2, 2, 4)
plot(tdx, dur(:, 4))
xlabel('time');
ylabel('VS');
axis('tight');
disp('press any key to continue')
pause
close all

disp('Perform principal components analysis')
disp(' ')
%perform singular value decomposition of the covariance matrix to determine
%the principal components.
[U, S, V] = svd(cov(dur));

in.fid = 1;
in.fmt = char('%12.4f');
tit = 'SD:';
mprintar(sqrt(S), in, tit);
tit = 'Loading Matrix:';
mprintar(U, in, tit);
disp('press any key to continue')
pause


disp('Estimate VAR(1)')
disp(' ')
nlag = 1;
res = var_est(dur, nlag);


disp(' ');
disp('***** Estimated VAR Model  *****');
disp(' ');
clear in
in.fid = 1;
in.fmt = char('%12.4f');
tit = 'AR';
strt = 1;
mprintar(res.phi(:, :, 2), in, tit, strt);
disp(' ')
tit = 'Constant';
mprintar(res.const', in, tit);

disp('Perform principal components analysis of the residuals')
disp(' ')
%perform singular value decomposition of the covariance matrix to determine
%the principal components.
[U, S, V] = svd(cov(res.resid));

in.fid = 1;
in.fmt = char('%12.4f');
tit = 'SD:';
mprintar(sqrt(S), in, tit);
tit = 'Loading Matrix:';
mprintar(U, in, tit);
disp('press any key to continue')
pause

h4p = [1, 0, -1, -1];
tit = 'h4p times phi(1):';
mprintar(h4p*res.phi(:, :, 2), in, tit);
disp('press any key to continue')
pause

tit = 'h4p times constant:';
mprintar(h4p*res.const, in, tit);
