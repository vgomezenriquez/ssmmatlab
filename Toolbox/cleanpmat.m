function [M, gM] = cleanpmat(M, tol)
%*************************************************************************
% This function cleans a polynomial matrix of small entries (in relation to its
% L1 norm) and if neccesary reduces its order afterwards
%
%  INPUTS:
%      M : (n x m x p) polynomial matrix
%    tol : tolerance value
%
% OUTPUTS:
%      M : (n x m x (gM+1)) polynomial matrix after elimination of small entries
%     gM : final order of the polynomial matrix
%*************************************************************************

if (nargin < 2)
    %     tol = 1.e-10;
    A = M(:);
    na = length(A);
    normaa = norm(A);
    tol = double(na) * normaa * eps(normaa);
end
Norm1 = max(abs(M(:)));
[n, m, p] = size(M);
tol = Norm1 * tol;
icual = (abs(M) < tol);
M(icual) = 0;
gM = 0;
for j = p:-1:2
    if (any(any(M(:, :, j))))
        gM = j - 1;
        M = M(:, :, 1:j);
        return
    end
end
M = M(:, :, 1);
return
