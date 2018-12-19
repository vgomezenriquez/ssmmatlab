%script file to simulate a VARMA model
%
% Let the VARMA model be
%
%   Phi(B) y_t = Theta(B) a_t,
%
% where B is the backshift operator, By_t = y_{t-1}.
%
clear 

l = 50;     %number of initial observations to be discarded in the
            %simulated series
N = 300;    %number of observations of the simulated series
s = 3;      %number of outputs
seed = 20;  %this is to always generate the same series

%polynomial matrices Phi and Theta
Phi(:, :, 1) = eye(s);
Theta(:, :, 1) = eye(s);
Theta(:, :, 2) = -[0.8, 0.4, 0.; -0.3, 0.6, 0.; 0., 0., -1.];

%covariance matrix of the a_t innovartions
S = [2.0, 0.5, 0.; 0.5, 1.0, 0.; 0., 0., .3];

%simulate v_t
[v, ferror] = varmasim(l, N, Phi, Theta, S, seed);

disp('cross correlation matrices')
%ccm matrices
lag = 6;
ic = 1;
str = mautcov(v, lag, ic);
disp('Correlation matrix at lag 0:')
disp(str.r0)
