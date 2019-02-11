%
%script file to illustrate two covariance factorization methods. The first
%one is Tunnicliffe-Wilson method and the second one uses the DARE.
%
clear

phi(:, :, 1) = eye(2);
th(:, :, 1) = eye(2);
th(:, :, 2) = [6.5, -2.; 15, -4.5];
Sigma = [4., 1.; 1., 1.];

nc = 2;
[c, ierror] = macgf(phi, th, Sigma, nc);
disp('Autocovariance matrices of lags 0 to 1:')
for i = 1:2
    disp(c(:, :, i))
end

disp('covariance factorization using Tunnicliffe-Wilson method')

[Omega, Theta, ierror, iter, normdif] = pmspectfac(c, 50)

disp('press any key to continue')
pause

disp('compute autocovariances of the invertible model')
phin(:, :, 1) = eye(2);
thn(:, :, 1) = eye(2);
thn(:, :, 2) = Theta;
Sigman = Omega;

nc = 2;
[cn, ierror] = macgf(phin, thn, Sigman, nc);
disp('Autocovariance matrices of lags 0 to 1:')
for i = 1:2
    disp(cn(:, :, i))
end
disp('press any key to continue')
pause

disp('covariance factorization using the DARE')
[Thetad, Omegad] = ssmspectfac(c)
