%script file for example3.7 in Tsay (2014)
%

phi(:, :, 1) = eye(2);
phi(:, :, 2) = -[0.816, -0.623; -1.116, 1.074];
phi(:, :, 3) = -[-0.643, 0.592; 0.615, -0.133];
th(:, :, 1) = eye(2);
th(:, :, 2) = -[0, -1.248; -0.801, 0];
Sigma = [4, 2; 2, 5]


nc = 3;
[c, ierror] = macgf(phi, th, Sigma, nc);
disp('Autocovariance matrices of lags 0 to 2:')
for i = 1:3
    disp(c(:, :, i))
end