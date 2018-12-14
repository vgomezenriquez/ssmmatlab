function [X, ierror] = sympmeq(phi, Lp)
%
% This function solves the symmetric polynomial matrix equation
%    X(z)phi'(z^{-1}) + phi(z)X'(z^{-1}) = Lp(z,z^{-1}),
% where Lp is a symmetric Laurent polynomial matrix of the form
%    Lp(z,z^{-1})= L'_pz^{-p}+.....+L'_1z^{-1}+L_0+L_1z+...L_pz^p,
% L_0 is a symmetric matrix, and p=degree(phi(z))= degree(X(z)).
%
% Input: phi = [n,n,p+1] matrix containing the phi polynomial matrix
%        Lp  = [n,n,p+1] matrix containing L_0+L_1z+...L_pz^p
% Output: X  = [n,n,p+1] matrix containing the solution X(z)
%        ierror = 0,1  a flag for errors in dimensions
%
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
%
[np, mp, kp] = size(phi);
[nl, ml, kl] = size(Lp);
ierror = 0;
if (np ~= nl | mp ~= ml | kp ~= kl)
    disp('wrong dimensions of phi or Lp in sympmeq');
    ierror = 1;
    return
end
s = np;
p = kp - 1;
X = zeros(size(phi));
n = floor(p/2);
ieyep = 0;
if ~any(any(phi(:, :, 1)-eye(np)))
    ieyep = 1;
end


A1 = [];
A2 = [];
for i = p:-1:n + 1
    A1 = [A1; vec(Lp(:, :, i+1))];
end
for i = n:-1:0
    A2 = [A2; vec(Lp(:, :, i+1))];
end
% A1
% A2
P1 = phi(:, :, 1);
Psi0 = eye(np); %compute the necessary Psi(z)=Phi^{-1}(z) weights
Psi = -phi(:, :, 2);
if (ieyep == 0), Psi0 = P1 \ Psi0;
    Psi = P1 \ Psi * Psi0;
end %first weight
for i = 2:p %rest of the weights
    A = [];
    for j = i:-1:2
        A = [A, phi(:, :, j)];
    end
    AA = -phi(:, :, i+1);
    if (ieyep == 0), AA = AA * Psi0;
    end
    Psix = AA - A * Psi;
    if (ieyep == 0), Psix = P1 \ Psix;
    end
    Psi = [Psi; Psix];
end
% Psi
psi = phi;
psi(:, :, 1) = Psi0;
for i = 2:p + 1
    psi(:, :, i) = Psi((i - 2)*s+1:(i - 1)*s, :);
end
% phi
% psi
% C = pmatmul(phi,psi)
%compute the blocks of T_p^{-1} and H_p
s2 = s * s;
nt11 = (p - n) * s2;
nt22 = (n + 1) * s2;
T11 = zeros(nt11, nt11);
T21 = zeros(nt22, nt11);
T22 = zeros(nt22, nt22);
H12 = zeros(nt11, nt22);
H21 = zeros(nt22, nt11);
H22 = zeros(nt22, nt22);

bphi = zeros(s2, s2, p+1);
fbphi = zeros(s2, s2, p+1);
Is2 = sparse(eye(s2));
Is = sparse(eye(s));
W = [];
for i = 1:s
    %  W=[W kron(Is,Is(:,i))];
    W = [W, kIv(Is(:, i))];
end

if (ieyep == 1)
    bphi(:, :, 1) = Is2;
    fbphi(:, :, 1) = W;
else
    bphi(:, :, 1) = kAI(psi(:, :, 1));
    %  bphi(:,:,1)=kron(psi(:,:,1),Is)
    %  fbphi(:,:,1)=kron(Is,phi(:,:,1))*W;
    fbphi(:, :, 1) = postmulW(kIA(phi(:, :, 1)), s);
end
for i = 2:p + 1
    bphi(:, :, i) = kAI(psi(:, :, i));
    %  bphi(:,:,i)=kron(psi(:,:,i),Is)
    %  fbphi(:,:,i)=kron(Is,phi(:,:,i))*W;
    fbphi(:, :, i) = postmulW(kIA(phi(:, :, i)), s);
end
% bphi
% fbphi

% p
% n
%compute T11, T21, T22, H12, H21 and H22
T = zeros((p - n)*s2);
H = T;
for i = 1:p - n
    for j = 1:i
        T((i - 1)*s2+1:i*s2, (j - 1)*s2+1:j*s2) = bphi(:, :, i-j+1);
        H((i - 1)*s2+1:i*s2, (p + 1 - j)*s2+1:(p + 2 - j)*s2) = fbphi(:, :, p+1-i+j);
    end
end
T11 = T(1:nt11, 1:nt11);
H12 = H(1:nt11, nt11+1:end);
clear T H
T = zeros((n + 1)*s2);
H = T;
for i = p - n + 1:p + 1
    for j = 1:i
        T((i - 1)*s2+1:i*s2, (j - 1)*s2+1:j*s2) = bphi(:, :, i-j+1);
        H((i - 1)*s2+1:i*s2, (p + 1 - j)*s2+1:(p + 2 - j)*s2) = fbphi(:, :, p+1-i+j);
    end
end
T21 = T(nt11+1:end, 1:nt11);
T22 = T(nt11+1:end, nt11+1:end);
H21 = H(nt11+1:end, 1:nt11);
H22 = H(nt11+1:end, nt11+1:end);
clear T H
% T11
% T21
% T22
% H12
% H21
% H22


%compute matrices to solve the linear system for X2
D12 = T11 * H12;
D21 = T22 * H21;
D22 = T21 * H12 + T22 * H22;
C1 = T11 * A1;
C2 = T21 * A1 + T22 * A2;
% D12
% D21
% D22
% C1
% C2

AX2 = eye(nt22) + D22 - D21 * D12;
BX2 = C2 - D21 * C1;
% AX2
% BX2
AX2 = sparse(AX2);
BX2 = sparse(BX2);
% if (ieyep == 1)
D = duplication(s); %introduce vech in X2
% D
ns = n * s2;
%  Ds=zeros(ns+s2,ns+s*(s+1)/2);
%  Ds(1:ns,1:ns)=eye(ns);
%  Ds(ns+1:end,ns+1:end)=D;
%  % Ds
%  % size(BX2)
%  % size(AX2)
%  % size(AX2*Ds)
%  Ds=sparse(Ds);
D = sparse(D);
%  X2=(AX2*Ds)\BX2;                %solve overdetermined system
AX2Ds = [AX2(:, 1:ns), AX2(:, ns+1:end) * D];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%the following line changed on 8-3-2010 to avoid singularity problems
%  X2=AX2Ds\BX2;                     %solve overdetermined system
X2 = pinv(full(AX2Ds)) * BX2;
%end of change
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  X2=Ds*X2;                       %restore vec in X2
X2 = [X2(1:ns); D * X2(ns+1:end)]; %restore vec in X2
% else
%  X2=AX2\BX2;
% end
X1 = C1 - D12 * X2;
Y = [X1; X2];
for i = 1:kp %pass from vec to matrix form
    for j = 1:np
        X(:, j, kp-i+1) = Y((i - 1)*s2+(j - 1)*s+1:(i - 1)*s2+j*s);
    end
end
% X
% % check that the solution is correct:
% disp('the following expression')
% Xp=pmmulbf(X,phi) + pmmulbf(phi,X)
% Lpp=zeros(np,np,2*p+1); Lpp(:,:,p+1)=Lp(:,:,1);
% for i=1:p
%  Lpp(:,:,i)=Lp(:,:,p+2-i)'; Lpp(:,:,p+1+i)=Lp(:,:,i+1);
% end
% disp('should be equal to')
% Lpp
% format long g
%  Xp-Lpp
% format short
