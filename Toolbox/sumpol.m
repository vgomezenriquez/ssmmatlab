function [C, ierrsumpol] = sumpol(A, B)
% This function sums two polynomials A and B, possibly of different
% degrees. The polynomials are assumed to be in reversed order, that is
% A(x) = A(1)*x^(p-1) + A(2)*x^(p-2) + ... + A(p-1)*x + A(p)

C = [];
ierrsumpol = 0;
[n, p] = size(A);
[n2, q] = size(B);
if (n ~= n2)
    ierrsumpol = 1;
    return
end
A = fliplr(A);
B = fliplr(B); %reverse order
maxpq = max([p, q]);
minpq = min([p, q]);
C = zeros(n, maxpq);
for k = 1:minpq
    C(:, k) = A(:, k) + B(:, k);
end
if (p > q)
    for k = minpq + 1:maxpq
        C(:, k) = A(:, k);
    end
else
    for k = minpq + 1:maxpq
        C(:, k) = B(:, k);
    end
end
C = fliplr(C); %reverse order