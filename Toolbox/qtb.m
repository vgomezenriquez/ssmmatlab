function QtB = qtb(A, B, p, q)
% given the output matrix A of function jqrt containing the
% J-unitary Housholder transformations and a matrix B, this
% function computes the product of Q'*B
%
%---------------------------------------------------
% USAGE: QtB=qtb(A,B,p,q)
% where:    A = an m x n matrix, m >= n
%           B = mb x nb matrix
%         p,q = integers such that J = diag(I_p,-I-q) is a signature
%               matrix with n=p+q
%---------------------------------------------------
% RETURNS: the product Q'*B
%---------------------------------------------------
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
[m, n] = size(A);
[mb, nb] = size(B);
if m ~= mb
    error('matrices A and B have not the same number of rows')
end
zero = 0.0d0;
QtB = zeros(m, nb); %this line changed 21-04-2008, (m,p)-->(m,nb)
for k = 1:nb
    b = B(:, k);
    for j = 1:n
        [pj, qj] = distnj(m, j, p, q);
        if (A(j, j) ~= zero)
            sum = jprod(A(j:end, j), b(j:end), pj, qj); % v'*x
            temp = sum / A(j, j); % v'*J*y*beta
            b(j:end) = b(j:end) - A(j:end, j) * temp; % x-(beta*v'*x)*v
        end
    end
    QtB(:, k) = b;
end
