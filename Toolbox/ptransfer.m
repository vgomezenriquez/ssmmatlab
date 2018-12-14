function [K, ierror] = ptransfer(A, B, n)
%
% This function computes the transfer function K(z)=A^{-1}(z)*B(z) up to
% the (n-1)-th term, that is, K(0),...,K(n-1)
% It is assumed that A(z) is square and that A(0) is nonsingular.
% Polynomial matrix B(z) can be nonsquare and, therefore, B(0) is not
% assumed to be the identity matrix.
%---------------------------------------------------
% USAGE: [K, ierror] = ptransfer ( A, B, n)
% where:    A = a square polynomial matrix with A(0) nonsingular
%           B = a polynomial matrix not necessarily square
%---------------------------------------------------
% RETURNS:
%           K = the first n weights of K(z)
%        ierror =1, dimension mismatch in A and B
%               =0, there are no errors on input
%---------------------------------------------------
%
% Copyright (c) 21 July 2003 by Victor Gomez
% Ministerio de Hacienda, Direccion Gral. de Presupuestos,
% Subdireccion Gral. de Analisis y P.E.,
% Alberto Alcocer 2, 1-P, D-34, 28046, Madrid, SPAIN.
% Phone : +34-915835439
% E-mail: VGomez@sepg.minhap.es
%
% The author assumes no responsibility for errors or damage resulting from the use
% of this code. Usage of this code in applications and/or alterations of it should
% be referenced. This code may be redistributed if nothing has been added or
% removed and no money is charged. Positive or negative feedback would be appreciated.
%

K = [];
ierror = 0;
[s1, s2, np] = size(A);
[s3, s4, nq] = size(B);
if (s1 ~= s2) | (s2 ~= s3)
    disp('wrong dimensions of A or B in ptransfer');
    ierror = 1;
    return
end
p = np - 1;
q = nq - 1;
K = repmat(0, [s1, s4, n]);
K0 = B(:, :, 1);
A0 = A(:, :, 1);
inv = 0;
if any(any(A0-eye(s1)))
    inv = 1;
    K0 = A0 \ K0;
end
K(:, :, 1) = K0;
for j = 1:n - 1
    if (q >= j)
        K(:, :, j+1) = B(:, :, j+1);
    end
    for i = 1:min(j, p)
        if (i < j)
            K(:, :, j+1) = K(:, :, j+1) - A(:, :, i+1) * K(:, :, j-i+1);
        else
            K(:, :, j+1) = K(:, :, j+1) - A(:, :, i+1) * K0;
        end
    end
    if inv == 1
        K(:, :, j+1) = A0 \ K(:, :, j+1);
    end
end
