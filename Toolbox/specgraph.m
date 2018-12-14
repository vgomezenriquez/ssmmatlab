function specgraph(th, phi, var)
%************************************************************************
% PURPOSE: This function graphs the spectrum of the ARIMA model
%
%          phi(B)y_t = th(B)*a_t
%
%---------------------------------------------------
% USAGE: specgraph(th,phi,var)
%
% Inputs: phi : a polynomial containing the AR part
%         th  : a polynomial containing the MA part
%        var  : a positive number, variance of the series model innovations
%
% Note: all of the previous polynomials are expressed in the same variable.
% The polynomials are given by an array like [ a_n, ... a_1, a_0], where
% the polynomial is a_0 + a_1*z + ... + a_n*z^n.
%
w = 0.001:0.01:pi;
z = exp(-sqrt(-1)*w);
%plot spectrum
sigma = sqrt(var);
h = polyval(th*sigma, z) ./ polyval(phi, z);
gainsq = abs(h).^2;
plot(w, gainsq);
axis([0, pi, 0, 3.])
