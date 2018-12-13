function p = cleanpol(p, tol)
% This function cleans a polynomial of small entries (in relation to its
% L1 norm) and if neccesary reduces its order afterwards
%  INPUTS:
%      p : (n x 1) polynomial coefficients
%    tol : tolerance value
%
% OUTPUTS:
%      p : polynomial after elimination of small entries

if (nargin < 2)
    na = length(p);
    normaa = norm(p);
    tol = double(na) * normaa * eps(normaa);
end
[n, m] = size(p);
icual = (abs(p) < tol); %returns the indices for which the condition holds
p(icual) = 0.;
p2 = 0;
%the following loop is to reduce the degree of the polynomial, if
%possible.
for j = 1:m
    if (icual(j))
        p2 = j;
    else
        break
    end
end
if p2 == m
    p = 0.;
else
    p = p(p2+1:end);
end
