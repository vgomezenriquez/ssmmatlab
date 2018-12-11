function [G, U, Dr, Nr, ierrglcd] = glcd(D, N, ivarmax)
% This function computes the Greatest Common Left Divisor G of a pair of
% polynomial matrices D and N of dimensions (n,n) and (n,m) respectively
% U is a unimodular matrix such that [D  N] U = [G 0] and G is lower triangular
%
% Moreover U = Vl -Nr  and if (D, N) was a left matrix fraction description of the
%              Ul  Dr      transfer function H then (Dr, Nr) is a right coprime
% matrix fraction description of H. If we compute inv(U), then
%  inv(U) = Dl  Nl and (Nl, Dl) is a left coprime matrix fraction description of H
%          -Ur  Vr

G = [];
U = [];
ierrglcd = 0;
if nargin < 3
    ivarmax = 0;
end
[n1, m1, p] = size(N);
[n, m, q] = size(D);
if (any([n, n]-[n1, m]))
    ierrglcd = 1; % dimension mismatch
    return
end
mpq = max(p, q);
NN = zeros(n1, m1, mpq);
DD = zeros(n, m, mpq);
NN(:, :, 1:p) = N;
DD(:, :, 1:q) = D;
M = cat(2, DD, NN);
% M = [D N];
% [G,U,rk,ierrLoTri] = LoTri(M, 1, iGaus, iHerm);
% We fill with the Identity matrix to guarantee full rank and triangularize
[G, U, ierrpmatri] = pmattrian(M, 1);
% M
% G
% U
% MU = pmatmul(M,U)
% pause
if (ierrpmatri > 0)
    ierrglcd = 2;
    return
end
G = G(1:n, 1:n, :);
Nr = -U(1:n, 1+n:n+m1, :);
Dr = U(n+1:n+m1, 1+n:n+m1, :);
if ivarmax == 1
    %check whether Dr(0) is the identity matrix. If not, change to VARMAX model
    A = Dr(:, :, 1);
    [nd, md, pd] = size(Dr);
    [nn, mn, pn] = size(Nr);
    if any(any(A-eye(m1)))
        for i = 2:pd
            Dr(:, :, i) = Dr(:, :, i) / A; %right division because it is a right MFD
        end
        for i = 1:pn
            Nr(:, :, i) = Nr(:, :, i) / A;
        end
        Dr(:, :, 1) = eye(m1);
    end
end
Nr = cleanpmat(Nr);
Dr = cleanpmat(Dr);
