function [Q, R, Indx, ierror] = housref(A, Q)
%
%---------------------------------------------------
% USAGE: [Q,R,Indx,ierror] = housref(A,Q)
% where:    A = an n x m matrix
%           Q = an n x n identity matrix
% Matrix A is reduced to row echelon form by means of Housholder
% transformations
%---------------------------------------------------
% RETURNS:
%           R = the upper triangular matrix in row echelon form such that A=Q*[R;0]
%           Q = a unitary matrix such that Q'*A = [R;0]
%        Indx = an index containing the l.i. and the l.d. columns in R
%        Indx(i)=0, i-th column is l.i.
%               =1, i-th column is l.d.
%        ierror =1, dimension mismatch in A and Q
%               =2,  Q is not the identity matrix on input
%               =0, there are no errors on input
%---------------------------------------------------
% If matrix Q is given on input, then the orthogonal Q matrix such that Q'*A = [R;0]
% is computed and given on output. If not, Q=[] on output.
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
ieyeq = 0;
ierror = 0;
R = [];
Indx = [];
if nargin < 2
    Q = [];
else
    [nq, mq] = size(Q);
    if (nq ~= m)
        disp('A and Q must have the same number of rows in housref')
        ierror = 1;
        return
    end
    if (sum(sum(Q-eye(m))) == 0)
        ieyeq = 1;
        Q1 = Q;
    else
        disp('Q must be the identity matrix on input in housref')
        ierror = 2;
        return
    end
end
mnmin = min(m, n);
if mnmin < 1
    return
end
% tol=double(m*n)*norm(A,'fro')*eps;
% normaa=norm(A); tol=double(m*n)*normaa*eps(normaa);
tol = max(size(A)) * eps(norm(A)); %Matlab's default

% We transpose the A matrix and apply the following function (Félix
% Aparicio, 2010).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The following function transforms an (n,m) matrix to lower triangular
% form using Householder reflections. A(input)*R =
% A(output) where R is an orthogonal matrix. See Stewart (1973),
% Introduction to Matrix Computations. Academic Press. pp. 230-235.
% A=A';  Pivots = zeros(m, 1); Indx=ones(1,n);
% for j=1:m
%  pj=0;
%  for k=j:n
%   if ( any( abs(A(k,j:m)) > tol) );
%    Pivots(j) = k;
%    pj = k; % We choose the first nonnull row of column j
%    break
%   end
%  end
%  if ( pj == 0)
%   R=A';
%   for i=1:m
%    if Pivots(i) > 0
%     Indx(Pivots(i))=0;
%    end
%   end
%   return %If all remaining rows of column j are null, we are done
%  end
%  x = A(pj,j:m); % this is the vector to transform to (ss,0,...,0)
%  xNorm = max(abs(x)); % xNorm is eta in 1) of page 234 of Stewart
%  v = x ./ xNorm; % v is v in 2) of page 234 of Stewart
%  if ( v(1) >= 0)
%   ss = sqrt(v*v');
%  else
%   ss = -sqrt(v*v');
%  end
%  u=v; u(1) = u(1)+ss;  % u is u in Stewart
%  pp=(u*u')/2;  % pp is pi
%  ss = xNorm*ss; % 6) in page 234 of Stewart
%  dim2 = m+1-j;
%  for k = pj:n
%   A(k,j:m) = A(k,j:m)-u.*((A(k,j:m)*u')/pp);
%  end
%  if ieyeq == 1
%   U = eye(dim2, dim2)-(1/pp).*(u'*u);
%   R2 = eye(m,m);  R2(j:m,j:m) = U;
%   Q = Q*R2;
%  end
% end
% R=A';
% for i=1:m
%  if Pivots(i) > 0
%   Indx(Pivots(i))=0;
%  end
% end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Original version, based on the function jqrt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R = zeros(m, n);
j = 1;
for i = 1:n
    if j <= m
        a = A(j:end, i);
        [na, ma] = size(a);
        Indx(i) = 0;
        [q, r] = jqrt(a, na, 0);
        Qta = qtb(q, A(j:end, i:end), na, 0);
        A(j:end, i:end) = Qta;
        if ieyeq == 1
            Qtq = qtb(q, eye(na), na, 0);
            Q2 = Q1;
            Q2(end-na+1:end, end-na+1:end) = Qtq;
            Q = Q2 * Q;
        end
        if jnorm(r, na, 0) <= tol
            Indx(i) = 1;
        else
            j = j + 1;
        end
    else
        Indx(i) = 1;
    end
end
if ieyeq == 1
    Q = Q';
end
R = A;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%