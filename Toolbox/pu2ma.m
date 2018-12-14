function [th, sigma2, ierrpu2ma] = pu2ma(p)
% This function obtains the moving average part corresponding
% to a finite covariance generating function that has been
% transformed into a polynomial in the variable U=z + z^(-1)

ierrpu2ma = 0;
rts = roots(p);
nrts = length(rts);
rth = [];
prod = 1.;
th = 1.;
for i = 1:nrts
    r = roots([1., -rts(i), 1.]);
    rth = [rth; r];
end
[aa, I] = sort(abs(rth), 'ascend');
rtho = rth(I);
n = length(I);
nn = n / 2;
for i = 1:nn
    prod = prod * (-rtho(i));
    th = conv(th, [-rtho(i), 1.]);
end
th = real(th);
prod = real(prod);
sigma2 = p(1) / prod;
