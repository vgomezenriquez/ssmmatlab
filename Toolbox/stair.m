function [a1, b1, k, u, C, v, ierror] = stair(a, b, tol, adj)
%
% [A1,B1,K,U,C,V,IERROR]=STAIR(A,B,TOL,ADJ) performs the staircase reduction of
% the pair (A,B). The transformed pair (A1,B1) has the typical staircase form
% (here with 4 "stairs"):
%
%                                     | X : * * * * * | } K(1)
%                                     |   : X * * * * | } K(2)
%         (B1:A1):=(U'*B*V:U'*A*U) =  |   :   X * * * | } K(3)
%                                     |   :     X * * | } K(4)
%                                     |   :         Z | }
%
% where each "stair" X has full row rank K(i) and has a lower triangular
% form (left adjusted) or an upper triangular form (right adjusted).
% In case adj is not given on input, the matrix V to adjust the X is not
% computed and the "stairs" X have just full row rank. The square matrix
% Z (if present) contains the uncontrollable modes of (A,B). The parameter
% tol is used to calculate the rank of the triangular matrices given by the
% qr algorithm.
%
% tol = [], the tolerence is computed by the program
%
%
% adj = 'l', lower triangular form (left adjusted)
%       'r', upper triangular form (right adjusted)
%
% This function is taken from VanDoren (2003) and modified by me on
% May 2008.
%
% Copyright (c) 2 January 2008 by Victor Gomez
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
[na, ma] = size(a);
[nb, mb] = size(b);
a1 = [];
b1 = [];
k = [];
u = [];
v = [];
ierror = 0;
C = []; % controllability indices
if (na ~= ma | na ~= nb)
    disp('wrong dimensions of a or b in stair');
    ierror = 1;
    return
end

if (nargin < 3)
    tol = [];
end

% otol=double(nb*nb)*max(norm(a,'fro'),norm(b,'fro'))*eps
if (isempty(tol))
    mnormab = max(norm(a), norm(b));
    tol = double(nb*nb) * mnormab * eps(mnormab);
end

n = nb;
m = mb; %k=zeros(1,n);
[q, r, e] = qr(b);
b1 = r * e';
a1 = q' * a * q;
u = q;
k(1) = rank(r);
k1 = 1;
k2 = k(1) + 1;
if (k2 <= n), b1(k2:n, 1:m) = zeros(n-k2+1, m);
end
for i = 2:n
    bh = a1(k2:n, k1:k2-1);
    ah = a1(k2:n, k2:n);
    [q, r, e] = qr(bh);
    q = [eye(k2-1), zeros(k2-1, n-k2+1); zeros(n-k2+1, k2-1), q];
    a1 = q' * a1 * q;
    a1(k2:n, k1:k2-1) = r * e';
    u = u * q;
    rk = rank(r, tol);
    if (k2 + rk <= n), a1(k2+rk:n, k1:k2-1) = zeros(n-k2-rk+1, k2-k1);
    end
    if rk == 0, break, end
    k(i) = rk;
    k1 = k2;
    k2 = k2 + rk;
    if k2 > n, break, end
end

%obtain controllability indices in descending order
C = zeros(1, m);
[nk, mk] = size(k);
for i = 1:mk
    for j = 1:k(i)
        C(j) = C(j) + 1;
    end
end
% C


if nargin == 4
    if (adj ~= 'l' & adj ~= 'r')
        disp('wrong argument adj in stair');
        ierror = 2;
        return
    end
    kmax = numel(k);
    nmax = sum(k);
    v = eye(mb);
    k3 = nmax;
    k2 = k3 - k(kmax) + 1;
    k1 = k2;
    if kmax > 1
        for i = kmax:-1:2
            k1 = k1 - k(i-1);
            if k2 > k1 + 1
                [q, r] = qr(pert(a1(k2:k3, k1:k2-1), adj));
                q = pert(q, adj)';
                a1(k1:k2-1, :) = q' * a1(k1:k2-1, :);
                a1(:, k1:k2-1) = a1(:, k1:k2-1) * q;
                u(:, k1:k2-1) = u(:, k1:k2-1) * q;
                a1(k2:k3, k1:k2-1) = pert(r, adj);
            end
            k3 = k2 - 1;
            k2 = k1;
        end
        if k3 > k2, b1(k2:k3, :) = q' * b1(k2:k3, :);
        end
    end
    [q, r] = qr(pert(b1(k2:k3, :), adj));
    q = pert(q, adj)';
    v = v * q;
    b1 = b1 * q;
    b1(k2:k3, :) = pert(r, adj); %upper triangular, right adjusted
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = pert(B, adj)
if adj == 'l'
    A = B'; %lower triangular, left adjusted
else
    A = fliplr(flipud(B')); %upper triangular, right adjusted
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
