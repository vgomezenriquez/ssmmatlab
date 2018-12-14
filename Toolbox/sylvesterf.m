function [Rd, fil, col] = sylvesterf(A, da, n, m, du);
% Given A, a polynomial matrix with n rows, m columns and degree da,
% this function computes its row-permuted Sylvester Matrix Rd. The permutation
% is done using the kronecker square commutation matrix of dimension n times
% da+du+1 Comm(n, da+du+1) (Actually this matrix is not computed and, equivalently
% but with fewer operations, the row permutation is done by the permat routine)

fil = (da + du + 1) * n;
col = (du + 1) * m;
Rd = zeros(fil, col);
filblo = (da + 1) * n;
bloque = zeros(filblo, m); %XXX Optimizar que bloque se calcule en LoTri una sola vez antes
for k = 1:da + 1 % de la primera llamada a esta rutina, pues no cambia de una a otra vez
    bloque(n*(k - 1)+1:n*k, :) = A(:, :, da+2-k);
end
for j = 1:du + 1
    Rd(n*(j - 1)+1:n*(j - 1)+filblo, m*(j - 1)+1:m*j) = bloque;
end
% P = commutation(n,da+du+1);
%Rd = P*Rd
Rd = permat(Rd, n, da+du+1);
