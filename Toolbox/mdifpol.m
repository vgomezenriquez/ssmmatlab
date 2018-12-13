function [Dr, Fd, ferror] = mdifpol(r, t, Pi1, Pid1)
%
%
%   This function obtains the differencing polynomial matrix correspondig
%   to the matrix polynomial pi(z) = eye(s) + pi_1*z + pi_2*z^2. Cases I(1)
%   and I(2) are considered.
%
% Input arguments:
%         r  : an integer, the rank of Pi1
%         t  : an integer, the rank of alphaor'*Pid(1)*betaor (Johansen)
%         Pi1: an s x s polynomial matrix, equal to -pi(1)
%        Pid1: an s x s polynomial matrix, equal to dpi(1)
% Output arguments:
%         Dr: a matrix polynomial, equal to the differencing matrix
%             polynomial
%         Fd: a matrix containing information for the parametrization of
%             the differencing matrix polynomial
%         ferror: a flag for errors
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

Dr = [];
Fd = [];
ferror = 0;

[s, ss] = size(Pi1);
if (s ~= ss)
    ferror = 1;
    disp('Pi1 should be square in mdifpol');
    return
end

if nargin < 4
    Pid1 = [];
end

J1 = t; %number of Jordan blocks of order one
J2 = s - r - t; %number of Jordan blocks of order two
prt = 0;
if (prt == 1)
    out = [s - r, r, J1, J2];
    fprintf(1, ' nr = %d, r = %d, J1 = %d, J2= %d\n', out);
end

twos = 2 * s;
threes = 3 * s;
[U, S, V] = svd(Pi1);
cont = 0;
M = [];
VJ = [];
I2 = 1;
I1 = 0;
if r > 0 %Pi1 is not zero
    %  Ssp5=S(1:r,1:r).^.5;
    %  alpha=U(:,1:r)*Ssp5; beta=V(:,1:r)*Ssp5;
    %  alphaor=U(:,r+1:end);
    betaor = V(:, r+1:end);
    %  aa= alphaor'*Pid1*betaor;
    Dr(:, :, 1) = eye(s);
    Dr(:, :, 2) = -betaor * pinv(betaor'*betaor) * betaor';
    Fd = betaor;
    return
    %  x0=betaor;  %eigenvectors corresponding to the eigenvalue one
    %  if J2 > 0
    %   if t == 0 %alphaor'*Pid(1)*betaor is zero
    %    A=alpha'*alpha*beta'; B=alpha'*Pid1*betaor;
    %    x1=A\B; %generalized eigenvectors, eigenvalue one
    %  % [x0(:,1:J2) x1] chains of order two, x0(:,J2+1:end) chains of order one
    %   else
    %    [Uo,So,Vo]=svd(alphaor'*Pid1*betaor);
    %    Ssp5=So(1:t,1:t).^.5;
    %    xi=Uo(:,1:t)*Ssp5; eta=Vo(:,1:t)*Ssp5;
    %    xietap=xi*eta'; etaor=Vo(:,t+1:end); x0etaor=x0*etaor;
    %    x0=[x0etaor x0*eta];
    %    x1=(alpha'*alpha*beta')\(alpha'*Pid1*x0etaor);
    %    % [x0etaor x1] chains of order two, x0*eta chains of order one
    %   end
    %   for i=1:J2
    %    cont=cont+1;
    %    v1=x0(:,cont); v2=x1(:,cont);
    %    M=[M v1 v2];
    %    S0v1=[v1; v1; v1]; S0v2=[v2; v2; v2];
    %    S1v1=[2.*v1; v1; zeros(size(v1))];
    %    VJ=[VJ S0v1 S0v2+S1v1];
    %   end
    %  end
    %  for i=1:J1
    %   cont=cont+1;
    %   v1=x0(:,cont);
    %   M=[M v1];
    %   S0v1=[v1; v1; v1];
    %   VJ=[VJ S0v1];
    %  end
else %Pi1 is zero
    if J2 > 0 %D(z)=(1-z)*(I-D_1*z), case I(2) reduced to I(1)
        [Uo, So, Vo] = svd(Pid1);
        x0 = Vo(:, t+1:end);
        for i = 1:J1
            cont = cont + 1;
            v1 = x0(:, cont);
            M = [M, v1];
            S0v1 = [v1; v1];
            VJ = [VJ, S0v1];
        end
        I2 = 0;
    else %D(z)=(1-z)^2*I, I(1) case
        Dr(:, :, 1) = eye(s);
        Dr(:, :, 2) = -eye(s);
        %   Fd=[ eye(s) zeros(s)];
        Fd = eye(s);
        %search is finished in I(1) case. We have obtained the maximum number of
        %unit roots (s).
        I1 = 1;
    end
end
%  M,VJ,I1,I2
%  [nv,mv]=size(VJ);
if I1 == 0 %case I(1) or I(2) is not finished yet
    [Q, R, Indx, ierror] = housref(VJ');
    [K, ierror] = nullref(R', Indx);
    [nk, mk] = size(K);
    Indk = zeros(1, nk);
    %determination of the leading entries in K. For each column, the leading
    %entry is stored in Indk.
    A = K(:);
    na = length(A);
    normaa = norm(A);
    tol = double(na) * normaa * eps(normaa);
    for i = 1:nk
        for j = mk:-1:1
            if (abs(K(i, j)) > tol)
                Indk(i) = j;
                break
            end
        end
    end
    %   K,Indk
    %determination of the minimum rows for the differencing polynomial
    indmax = 0;
    col = 0;
    for i = 1:s
        minc = i + threes;
        indm = 0;
        for j = 1:nk
            indkj = Indk(j);
            if (mod(indkj-i, s) == 0) && (indkj < minc) && (indkj > col)
                minc = indkj;
                indm = j;
                col = indkj;
                if minc > indmax
                    indmax = minc;
                end
            end
        end
        MD(i, :) = K(indm, :);
    end
    %   MD, indmax
    Dr(:, :, 1) = eye(s);
    if (I2 == 0) %D(z)=(1-z)*(I-D_1*z), case I(2) reduced to I(1)
        D1 = MD(:, 1:s);
        Dr(:, :, 2) = D1;
        Da(:, :, 1) = eye(s);
        Da(:, :, 2) = -eye(s);
        Dr = pmatmul(Da, Dr);
        Fd = -[Dr(:, :, 2), Dr(:, :, 3)];
    else %the differencing matrix polynomial is MD= [D_2 D_1 I]
        if (indmax == threes)
            D1 = MD(:, s+1:twos);
            D2 = MD(:, 1:s);
            Dr(:, :, 2) = D1;
            Dr(:, :, 3) = D2;
            Fd = -[D1, D2];
        else %the differencing matrix polynomial is MD= [D_1 I] (D2=0)
            D1 = MD(:, 1:s); %eig(D1), rank(eye(s)+D1)
            Dr(:, :, 2) = D1;
            Fd = [-D1, zeros(s)];
        end
    end
end
