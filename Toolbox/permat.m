function [B] = permat(A, n, m);
% This function permutes the rows of a matrix A (that has nm rows) as if we
% premultiplied A by Knm, the permutation matrix of parameters n and m
%
%permutation matrix: P*vec(A) = vec(A')
%

B = [];
[fil, col] = size(A);
% if ( fil-n*m)
%     ierrCommute = 1;
%     return
% end
B = zeros(fil, col);
icual = 1;
for i = 1:n
    for j = 1:m
        B(icual, :) = A((j - 1)*n+i, :);
        icual = icual + 1;
    end
end
