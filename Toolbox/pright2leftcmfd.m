function [phie, thetae, kro, ierror] = pright2leftcmfd(phir, thetar, np, kro)
%
% This function computes a left coprime MFD given a right MFD. That is,
% given thetar(z)*phir^{-1}(z), a left coprime MFD  is computed such that
% phie^{-1}(z)*thetae(z) = thetar(z)*phir^{-1}(z).
% It is assumed that phir(z) is square and that phir(0) is nonsingular.
% Polynomial matrix thetar(z) can be nonsquare and, therefore, thetar(0) is not
% assumed to be the identity matrix.
%---------------------------------------------------
% USAGE: [phie,thetae,kro,ierror] = pright2leftcmfd(phir,thetar,np,kro)
% where:    phir   = a k x k polynomial matrix with phir(0) nonsingular
%           thetar = a k x m polynomial matrix
%           np     = un upper bound for the Kronecker indices
%           kro    = a 1 x k vector containing the Kronecker indices
%---------------------------------------------------
% RETURNS:
%           phie  = the AR echelon polynomial matrix
%          thetae = the MA echelon polynomial matrix
%           kro   = a 1 x k vector containing the Kronecker indices
%        ierror =1, dimension mismatch in phir and thetar
%               =0, there are no errors on input
%---------------------------------------------------
% If kro is not input, the function uses functions housref and nullref on
% the augmented Sylvester matrices constructed with phir and thetar to
% compute the Kronecker indices. If kro is input, a system of linear
% equations based on an appropriate augmented Sylvester matrix is solved.
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

phie = [];
thetae = [];
ierror = 0;
[p1, p2, p3] = size(phir);
[t1, t2, t3] = size(thetar);
if (p1 ~= p2) | (p1 ~= t2)
    disp('wrong dimensions of phi or theta in pright2leftcmfd');
    ierror = 1;
    return
end
if nargin < 3
    kro = [];
    ierror = 2;
    return
end
if nargin < 4
    ikro = 0;
else
    ikro = 1;
end
s = t1;
m = p2; %t1=s, p1=m=p2;
nq = np;

d = max(p3, t3);
bphir = zeros(p1, p2, d);
bthetar = zeros(t1, t2, d);
bphir(:, :, 1:p3) = phir;
bthetar(:, :, 1:t3) = thetar;
% check that indeed D(0)=I. If not, change D and N.
A = bphir(:, :, 1);
if any(any(A-eye(p1)))
    bphir(:, :, 1) = eye(p1);
    for i = 2:p3
        bphir(:, :, i) = bphir(:, :, i) / A; %right division because it is a right MFD
    end
    for i = 1:t3
        bthetar(:, :, i) = bthetar(:, :, i) / A;
    end
end

D = [];
N = [];
for i = d:-1:1
    D = [D, bphir(:, :, i)];
    N = [N, bthetar(:, :, i)];
end
% d,D,N

Indn = 1:s;
spm = s + m;
N0 = N;
DN = D; % elements of augmented Sylvester matrix

if ikro == 0 % Kronecker indices are unknown
    iser = 1;
    cont = 0;
    ckro = 0;
    kro = zeros(1, s);
    Ili = 1:m; % index for independent rows
    P = zeros(s);
    T = zeros(s, m);
    phie(:, :, 1) = P;
    thetae(:, :, 1) = T;
    while iser
        cont = cont + 1;
        [n0, m0] = size(N0);
        N00 = N0;
        Indn0 = Indn;
        Indn = [];
        N0 = [];
        cspm = cont * spm;
        cspmms = cspm - s;
        for i = 1:n0 % search last N block from top to bottom
            DNi = [DN; N00(i, :)]; % augment DN matrix with row of last N block
            [Q, R, Indx, ierror] = housref(DNi'); % test whether last row is l.d.
            if Indx(end) == 1 % l.d. row found in N block
                ckro = ckro + 1;
                kro(Indn0(i)) = cont - 1; % adjust Kronecker indices
                Idelx = Indn0(i); % create index for echelon row
                [K, ierror] = nullref(R', Indx); % solve by back substitution to find coefficients
                KK = zeros(1, cspm); % prepare vector with zeros for echelon row
                lastc = cspmms + Idelx; % position of the one in phi(0) in echelon row
                ix = [Ili, lastc]; % index for the columns in the row that are not zero
                %the following line modified on February, 3rd, 2009
                KK(:, ix) = K(end, :); % row for echelon form
                for j = 1:cont % fill in phi and theta rows
                    tpr = KK(:, (j - 1)*spm+1:j*spm);
                    phie(Idelx, :, cont-j+1) = tpr(m+1:end);
                    thetae(Idelx, :, cont-j+1) = -tpr(1:m);
                end
            else
                DN = [DN; N00(i, :)]; % add row in DN matrix
                Ili = [Ili, cspmms + Indn0(i)]; % adjust index for independent rows
                Indn = [Indn, Indn0(i)]; % store remaining l.i. indices only
                N0 = [N0; N00(i, :)]; % store l.i. rows of last N block only
            end
        end
        if (cont == np) | (ckro == s) % stop criteria
            iser = 0;
        end
        [nd, md] = size(DN); % augment Sylvester matrix with l.i. rows in N block
        DN = [[DN, zeros(nd, m)]; [zeros(m, cont*m), D]];
        [n0, m0] = size(N0);
        [ndn, mdn] = size(DN); % N0 contains the remaining l.i. rows in last N block
        N0 = [zeros(n0, mdn-m0), N0]; % adjust size of N0
        cspm = cont * spm;
        Ili = [Ili, cspm + 1:cspm + m]; % adjust index for independent rows
    end
else % Kronecker indices are known
    iser = 1;
    cont = 0;
    ckro = 0;
    Ili = 1:m; % index for independent rows
    P = zeros(s);
    T = zeros(s, m);
    phie(:, :, 1) = P;
    thetae(:, :, 1) = T;
    ldrow = zeros(size(kro));
    for i = 1:s % find row numbers for l.d. rows
        ldrow(i) = kro(i) * spm + m + i;
    end
    %  ldrow
    while iser
        cont = cont + 1;
        [n0, m0] = size(N0);
        N00 = N0;
        Indn0 = Indn;
        Indn = [];
        N0 = [];
        cspm = cont * spm;
        cspmms = cspm - s;
        for i = 1:n0 % search last N block from top to bottom
            DNi = [DN; N00(i, :)]; % augment DN matrix with row of last N block
            ipld = cspmms + Indn0(i); % index for added row
            if ldrow(Indn0(i)) == ipld % added row is l.d.
                [Q, R, Indx, ierror] = housref(DNi');
                ckro = ckro + 1; % adjust Kronecker index count
                Idelx = Indn0(i); % create index for echelon row
                [K, ierror] = nullref(R', Indx); % solve by back substitution to find coefficients
                KK = zeros(1, cspm); % prepare vector with zeros for echelon row
                lastc = cspmms + Idelx; % position of the one in phi(0) in echelon row
                ix = [Ili, lastc]; % index for the columns in the row that are not zero
                KK(:, ix) = K; % row for echelon form
                for j = 1:cont % fill in phi and theta rows
                    tpr = KK(:, (j - 1)*spm+1:j*spm);
                    phie(Idelx, :, cont-j+1) = tpr(m+1:end);
                    thetae(Idelx, :, cont-j+1) = -tpr(1:m);
                end
            else
                DN = [DN; N00(i, :)]; % add row in DN matrix
                Ili = [Ili, ipld]; % adjust index for independent rows
                Indn = [Indn, Indn0(i)]; % store remaining l.i. indices only
                N0 = [N0; N00(i, :)]; % store l.i. rows of last N block only
            end
        end
        if (cont == np) | (ckro == s) % stop criteria
            iser = 0;
        end
        [nd, md] = size(DN); % augment Sylvester matrix with l.i. rows in N block
        DN = [[DN, zeros(nd, m)]; [zeros(m, cont*m), D]];
        [n0, m0] = size(N0);
        [ndn, mdn] = size(DN); % N0 contains the remaining l.i. rows in last N block
        N0 = [zeros(n0, mdn-m0), N0]; % adjust size of N0
        Ili = [Ili, cspm + 1:cspm + m]; % adjust index for independent rows
    end
end
if ckro ~= s
    ierror = 3; %no solution found, probably we have to increase np.
end
