function [Q, R] = jqrt(A, p, q)
%
%---------------------------------------------------
% USAGE: [Q,R] = jqrt(A,p,q)
% where:    A = an n x m matrix with n >=m
%           p,q = integers such that J=diag(I_p,-I_q) is a signature matrix
% It is assumed that A'*J*A = R'*J*R, where R is an upper triangular
% matrix. Then, there exists a J-unitary matrix Q such that A=Q*[R;0]
% To construct the matrix Q, J-unitary Housholder transformations are used
%---------------------------------------------------
% RETURNS:
%           R = the upper triangular matrix such that A=Q*[R;0]
%           Q = a factored form of a J-unitary matrix such that QJQ'=Q'JQ=J and
%               Q'*A = [R;0]
%           If the product Q'*B is desired for some matrix B, call function
%           qtb of this library (QtB=qtb(A,B,p,q)).
%---------------------------------------------------
%        A is an m by n array. on input A contains the matrix for
%          which the qr factorization is to be computed. on output
%          the strict upper trapezoidal part of A contains the strict
%          upper trapezoidal part of R, and the lower trapezoidal
%          part of a contains a factored form of Q.
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
[m, n] = size(A);
if m ~= p + q
    error('number of rows of A nonequal to p + q in jqrt')
end
one = 1.0d0;
p05 = 5.0d-2;
zero = 0.0d0;
%
%  epsmch is the machine precision.
%
%       epsmch = 2.225073858507201d-14;
%       epsmch=double(m*n)*norm(A,'fro')*eps;  %changed on 10/06/2008
epsmch = eps;
%
%  compute the initial column norms and initialize several arrays.
%
rdiag = zeros(n, 1);
wa = zeros(n, 1);
for j = 1:n
    ajnorm = jnorm(A(:, j), p, q);
    if ajnorm < zero
        error('J-norm negative in jqrt')
    end
    rdiag(j) = ajnorm;
    wa(j) = ajnorm;
end
%
%  reduce a to r with householder transformations.
%
minmn = min(m, n);
for j = 1:minmn
    %
    %  compute the householder transformation to reduce the
    %  j-th column of A to a multiple of the j-th unit vector.
    %
    [pj, qj] = distnj(m, j, p, q);
    ajnorm = jnorm(A(j:end, j), pj, qj);
    if ajnorm < zero
        error('J-norm negative in jqrt')
    end
    if (ajnorm ~= zero)
        if (A(j, j) < zero)
            ajnorm = -ajnorm;
        end
        A(j:end, j) = A(j:end, j) / ajnorm;
        A(j, j) = A(j, j) + one; % the transformation is Hx=(I-beta*v*v'*J)x, where v=x/s+e_1,
        % s=sqrt(x'*J*x) and beta=1/v_1.
        %
        % apply the transformation to the remaining columns
        % and update the norms.
        %
        jp1 = j + 1;
        if (n >= jp1)
            for k = jp1:n
                sum = jprod(A(j:end, j), A(j:end, k), pj, qj); % v'*J*y
                temp = sum / A(j, j); % (v'*J*y)*beta
                A(j:end, k) = A(j:end, k) - temp * A(j:end, j); % y-(v'*J*y*beta)*v
                if (rdiag(k) ~= zero)
                    temp = A(j, k) / rdiag(k);
                    rdiag(k) = rdiag(k) * sqrt(max(zero, one-temp^2));
                    if (p05 * (rdiag(k) / wa(k))^2 <= epsmch)
                        %the following line changed on 17-3-2010.
                        [pjj, qjj] = distnj(m, jp1, p, q);
                        rdiag(k) = jnorm(A(jp1:end, k), pjj, qjj);
                        wa(k) = rdiag(k);
                    end
                end
            end
        end
    end
    rdiag(j) = -ajnorm;
end

% compute matrix R
R = A;
for i = 1:n
    R(i, i) = rdiag(i);
end
for i = 2:m
    k = min(n, i-1);
    R(i, 1:k) = zeros(1, k);
end

Q = A; % matrix Q is stored in A as a sequence of Housholder transformations


% % compute matrix Q
% Q=eye(m);
% for k=1:m
%  b=Q(:,k);
%  for j=1:n
%   [pj,qj]=distnj(m,j,p,q);
%   if (A(j,j) ~= zero)
%    sum=jprod(A(j:end,j),b(j:end),pj,qj); % v'*x
%    temp=-sum/A(j,j);                     % -beta*v'*x
%    b(j:end)=b(j:end)+A(j:end,j)*temp;    %x-(beta*v'*x)*v
%   end
%  end
%  Q(:,k)=b;
% end
% Q=Q';
