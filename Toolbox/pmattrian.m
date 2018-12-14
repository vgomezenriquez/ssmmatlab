function [T, U, ierrpmatri] = pmattrian(A, full, iHerm)
% This function computes T, a lower triangular form of a full column rank polinomial
% matrix A. U is a unimodular matrix such that AU = T and rk = rank(A) can be used on
% output to check if A really had full rank. If it hadn´t, on output T and U are empty.
% See Henrion and Sebek (1999) Reliable Numerical Methods for Polynomial Matrix
% Triangularization, IEEE Trans. Aut. Control 44-3 pp. 497-508
% Author Felix Aparicio-Perez, Instituto Nacional de Estadistica, Spain

[n, m, p] = size(A);
T = [];
U = [];
ierrpmatri = 0; %Initialization

if (n == 1) & (m == 1)
    T = A;
    U = 1.;
    return
end

M = A(:);
nm = length(M);
normam = norm(M);
tol = double(nm) * normam * eps(normam);

% [n,m,p] = size(A);
% T=[]; U=[]; rk=0; ierrpmatri = 0; %Initialization

% If full ~= 0 we augment A with the (m,m) identity matrix in order to ensure full
% column rank of the resulting matrix; if full = 0 no completion is made and we
% assume the risk of A not having full column rank
if (full)
    A = cat(1, A, zeros(m, m, p));
    A(n+1:n+m, :, 1) = eye(m, m);
    n0 = n;
    n = n + m;
end
dr = zeros(n, 1);
ds = zeros(m, 1); %Initialization of degrees of rows and columns
% Compute the degree of A, an upper bound dumax of du and the row and
% column degrees
for i = 1:n
    for k = p:-1:1
        if (any(abs(A(i, :, k)) > tol));
            dr(i) = k - 1;
            break
        end
    end
end
da = max(dr);
for j = 1:m
    for k = p:-1:1
        if (any(abs(A(:, j, k)) > tol));
            ds(j) = k - 1;
            break
        end
    end
end
if (da == -1)
    ierrpmatri = 1; %In this case matrix A is null
    return
end
sdr = sort(dr, 'descend');
sds = sort(ds, 'descend');
dumax = min([sum(sdr(1:m-1)), sum(sds(1:m-1))]);
iGaus = 0;
iHerm = 0;
% Iterate from du=0 to dumax ( best if by bipartition, improvement XXX)
for du = 0:dumax
    dtot = da + du;
    %    du
    [Rd, fil, col] = sylvesterf(A, da, n, m, du);
    %   [Um,R] = qr(Rd');
    if (full == 0) %We didn´'t fill with Im, so we need to have Um computed
        %             [Tm, ierrHousePar, Um] = HousePar(Rd);
        [Um, Tmt, Indx, ierror] = housref(Rd', eye(col));
        Tm = Tmt';
    else % We filled with Im, so we'll get U from the lower part of T
        Um = [];
        %            [Tm, ierrHousePar] = HousePar(Rd);
        [Um, Tmt, Indx, ierror] = housref(Rd');
        Tm = Tmt';
    end
    [SigBar] = mshape(Tm); %We compute the shape of Tm
    %We compute the triangular shape, if it is defined
    [Sig, ColsUT, Cs, Ncs] = csigsets(SigBar, dtot, du, n, m);
    if (~isempty(Sig)) %We have found a triangular shape
        % We compute the Unimodular polynomial matrix U
        % and the polynomial matrix T from Um and Tm
        [U, T] = copmut(Um, Tm, ColsUT, da, du, n, m, full);
        % We clean unnecessary zeros from U and T
        U = cleanpmat(U);
        T = cleanpmat(T);
        return
    end
end
% Something strange happened and we didnn't find a triangular shape
ierrpmatri = 2;
