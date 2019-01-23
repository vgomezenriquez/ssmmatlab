function [KKP, PT, hd, Md] = smoothgen(y, X, Z, G, W, T, H, ins, i, mucd, U, C, D)
%
%
%        This function applies the augmented Kalman filter and smoother
%        to the series y corresponding to the model
%
%        y_t = X_t*beta + Z_t*alpha_t + G_t*epsilon_t
%        alpha_{t+1}= W_t*beta + T_t*alpha_t + H_t*epsilon_t
%
%        where epsilon_t is (0,sigma^2I),
%
%        with initial state
%
%        alpha_1= c + W_0*beta + a_1 + A_1*delta
%
%        where c is (0,Omega) and delta is (0,kI) (diffuse). A single
%        collapse is applied to get rid of the diffuse component.
%
%        It is desired to smooth the vector
%
%         Y_t = U_t*\beta + C_t*alpha_t + D_t*epsilon_t
%
%        Input parameters:
%        y:     an (n x p) matrix of observations;
%        X    : an (n*p x nbeta) matrix containing the X_t matrices;
%               a  (p x nbeta) if it is time invariant;
%               it can be []
%        Z    : an (n*p x nalpha) matrix containing the Z_t matrices;
%               a  (p x nalpha) matrix if it is time invariant
%        G    : an (n*p x nepsilon) matrix containing the G_t matrices;
%               a  (p x nepsilon) matrix if it is time invariant
%        W    : an (n*nalpha x nbeta) matrix containing the W_t matrices;
%               an (nalpha x nbeta) matrix if it is time invariant;
%               it can be []
%        T    : an (n*nalpha x nalpha) matrix containing the T_t matrices;
%               an (nalpha x nalpha) matrix if it time invariant
%        H    : an (n*nalpha x nepsilon) matrix containing the H_t matrices;
%               an (nalpha x nepsilon) if it is time invariant
%        ins: an nalpha x (cc+cw0+ca1+cca1) matrix containing the initial
%             state information, according to array i below
%        i    : a  1 x 4 array containing 4 integers, i=[cc cw0 ca1 cca1],
%             where
%             cc   = nalpha if c is not missing (0 if c missing)
%             cw0  = number of columns in W_0 (0 if W_0 missing)
%             ca1  = 1 if a_1 is not missing (0 if a_1 missing)
%             cca1 = number of columns in A_1 (0 if A_1 missing)
%       mucd: an integer, the dimension of Y_t
%        U  : an (n*mucd x nbeta) matrix containing the U_t matrices;
%             an (mucd x nbeta) if it is time invariant
%             it can be []
%        C  : an (n*mucd x nalpha) matrix containing the C_t matrices;
%             an (mucd x nalpha) if it is time invariant
%        D  : an (n*mucd x nepsilon) matrix containing the D_t matrices;
%             an (mucd x nepsilon) if it is time invariant
%
%        Output parameters:
%        KKP : an (n x mucd) matrix containing the estimated Y_{t|n}
%        PT  : an (mucd*n x mucd) matrix containing the Mse of Y_{t|n}
%        hd  : the (delta',beta')' estimate
%        Md  : the Mse of hd
%
% Copyright (c) 02 December 2005 by Victor Gomez
% Ministerio de Economia y Hacienda, Direccion Gral. de Presupuestos,
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
%
% system dimensions
%
[n, p] = size(y);
[mx, nx] = size(X);
[mz, nalpha] = size(Z);
[mg, neps] = size(G);
[mw, nw] = size(W);
[mt, nt] = size(T);
[mh, nh] = size(H);
[mu, nu] = size(U);
[mmc, nnc] = size(C);
[md, nd] = size(D);
[mi, ni] = size(ins);
[mc, nc] = size(i);
nbeta = max(nx, nw);
cc = i(1);
cw0 = i(2);
ca1 = i(3);
cca1 = i(4);
%
% check for inconsistencies
%
np = n * p;
nnalpha = n * nalpha;
nmucd = n * mucd;
if nbeta > np
    disp('too many regressors in smoothgen')
    return
end
if mi ~= nalpha
    disp('the number of rows in ins is incorrect')
    return
end
if (mmc ~= md)
    disp('the matrices C and D should have the same number of rows')
    return
end
if (mu > 0) && (mu < md)
    disp('the matrix U hass less rows than the matrix D')
    return
end
if (mu > 0) && (mu < mmc)
    disp('the matrix U has less rows than the matrix C')
    return
end
if mx > p || mz > p || mg > p || mw > nalpha || mt > nalpha ...
        || mh > nalpha || mu > mucd || mmc > mucd || md > mucd
    %
    % system matrices are time varying
    %
    if mx > p && mx ~= np
        disp('the number of rows in X is incorrect')
        return
    end
    if mz > p && mz ~= np
        disp('the number of rows in Z is incorrect')
        return
    end
    if mg > p && mg ~= np
        disp('the number of rows in G is incorrect')
        return
    end
    if mw > nalpha && mw ~= nnalpha
        disp('the number of rows in W is incorrect')
        return
    end
    if mt > nalpha && mt ~= nnalpha
        disp('the number of rows in T is incorrect')
        return
    end
    if mh > nalpha && mh ~= nnalpha
        disp('the number of rows in H is incorrect')
        return
    end
    if mu > mucd && mu ~= nmucd
        disp('the number of rows in U is incorrect')
        return
    end
    if mmc > mucd && mmc ~= nmucd
        disp('the number of rows in C is incorrect')
        return
    end
    if md > mucd && md ~= nmucd
        disp('the number of rows in D is incorrect')
        return
    end
end
if nw ~= nx && nw ~= 0 && nx ~= 0
    disp('the number of columns in W and X should agree')
    return
end
if nw ~= nu && nw ~= 0 && nu ~= 0
    disp('the number of columns in W and U should agree')
    return
end
if nt ~= nalpha
    disp('the number of columns in T and Z should agree')
    return
end
if nnc ~= nalpha
    disp('the number of columns in T and C should agree')
    return
end
if nh ~= neps
    disp('the number of columns in H and G should agree')
    return
end
if nd ~= neps
    disp('the number of columns in D and G should agree')
    return
end
if mi ~= nalpha
    disp('the number of rows in ins and the number of columns in')
    disp(' Z should agree')
    return
end
if nc ~= 4 || mc ~= 1
    disp('c should be a 1 x 4 matrix')
    return
end
if nbeta == 0 && cw0 ~= 0
    disp('the number of rows in W_0 should be zero')
    return
end
if ni ~= cc + cw0 + ca1 + cca1
    disp('the number of columns in ins should be ')
    disp(cc+cw0+ca1+cca1)
    return
end
%
% initial covariance matrix
%
if cc == 0
    P = zeros(nalpha, nalpha);
else
    P = ins(:, 1:cc);
end
%
% initial state
%
if cw0 == 0
    if nbeta == 0
        W0 = [];
    else
        W0 = zeros(nalpha, nbeta);
    end
else
    W0 = ins(:, cc+1:cc+cw0);
end
%W0
if ca1 == 0
    a1 = zeros(nalpha, 1);
else
    a1 = ins(:, cc+cw0+1:cc+cw0+ca1);
end
if cca1 == 0
    ndelta = 0;
    aa1 = [];
else
    ndelta = cca1;
    aa1 = ins(:, cc+cw0+ca1+1:cc+cw0+ca1+cca1);
end
A = [-aa1, -full(W0), a1];
%
% augmented Kalman filter
%
nd = ndelta;
ndb = ndelta + nbeta;
ndb1 = ndb + 1;
nb1 = nbeta + 1;
AV = zeros(n*p, ndb1);
AHD = zeros(nd, ndb1);
FST = zeros(n*p, p);
nkkp = max(nalpha, mucd);
KKP = zeros(n*p, nkkp);
AAT = sparse((n + 1)*nalpha, ndb1);
PT = zeros((n + 1)*nkkp, nkkp);
SQT = zeros(ndb+p, ndb1);
AAT(1:nalpha, :) = A;
PT(1:nalpha, 1:nalpha) = P;
iti = 0;
collps = 0;
if nd == 0
    collps = 1;
end
ifg = 0;
ifh = 0;
if mg <= p
    GG2i = G * G';
    ifg = 1;
end
if mh <= nalpha
    HH2 = H * H';
    if ifg == 1
        HGi = H * G';
    end
    ifh = 1;
end
%
for i = 1:n
    %   i
    ip = (i - 1) * p + 1:i * p;
    ia = (i - 1) * nalpha + 1:i * nalpha;
    if mx > p
        XX = X(ip, :);
    else
        XX = X;
        if nbeta > 0 && isempty(X) == 1
            XX = zeros(p, nbeta);
        end
    end
    if mz > p
        ZZ = Z(ip, :);
    else
        ZZ = Z;
    end
    if mg > p
        GG = G(ip, :);
    else
        GG = G;
    end
    if mw > nalpha
        WW = W(ia, :);
    else
        WW = W;
        if nbeta > 0 && isempty(W) == 1
            WW = zeros(nalpha, nbeta);
        end
    end
    if mt > nalpha
        TT = T(ia, :);
    else
        TT = T;
    end
    if mh > nalpha
        HH = H(ia, :);
        HH2 = HH * HH';
    else
        HH = H;
    end
    YY = y(i, :)';
    %
    % check for missing values
    %
    nn = isnan(YY);
    miss = sum(nn);
    if miss > 0 && miss < p
        idn = find(~nn);
        if ~isempty(XX)
            XX = XX(idn, :);
        end
        YY = YY(idn);
        if ~isempty(GG)
            GG = GG(idn, :);
            GG2 = GG * GG';
            HG = HH * GG';
        end
        if ~isempty(ZZ)
            ZZ = ZZ(idn, :);
        end
        %
        %tskf
        %
        V = [zeros(p-miss, ndelta), full(XX), YY] - ZZ * A;
        F = ZZ * P * ZZ' + GG2;
        [CF, cp] = chol(F);
        if cp == 0
            Frt = CF';
            DD = Frt \ V;
            K = (TT * P * ZZ' + HG) / F;
        elseif cp == 1 %F is zero
            DD = zeros(p-miss, ndb1);
            K = zeros(nalpha, p-miss);
            V = zeros(size(V));
            miss = p;
        else
            error('innovations covariance singular in smoothgen')
        end
    elseif miss == p
        DD = zeros(p, ndb1);
        K = zeros(nalpha, p);
        V = zeros(p,size(A,2));
    else
        if ifg == 0
            GG2 = GG * GG';
            HG = HH * GG';
        else
            GG2 = GG2i';
            if ifh == 1
                HG = HGi;
            else
                HG = HH * GG';
            end
        end
        V = [zeros(p, ndelta), full(XX), YY] - ZZ * A;
        F = ZZ * P * ZZ' + GG2;   
        [CF, cp] = chol(F);
        if cp == 0
            Frt = CF';
            DD = Frt \ V;
            K = (TT * P * ZZ' + HG) / F;
        elseif cp == 1 %F is zero
            DD = zeros(p, ndb1);
            K = zeros(nalpha, p);
            V = zeros(size(V));
            miss = p;
        else
            error('innovations covariance singular in smoothgen')
        end
    end
    %
    % SQT updating
    %
    if ndb > 0
        [~, SQT] = qr([DD; SQT(1:ndb, :)]);
        %     else
        %      SQT=DD;
    end
    if miss == p
        A = [zeros(nalpha, ndelta), -full(WW), zeros(nalpha, 1)] + TT * A;
        P = TT * P * TT' + HH2;
        AV(ip, :) = [zeros(p, nd-ndelta), zeros(p, ndb1)];
        FST(ip, :) = zeros(p, p);
        KKP(ip, 1:nalpha) = zeros(p, nalpha);
    else
        A = [zeros(nalpha, ndelta), -full(WW), zeros(nalpha, 1)] + TT * A + K * V;
        P = TT * P * (TT - K * ZZ)' + HH2 - HG * K';
        if miss == 0
            idn = 1:p;
        end
        DDm = zeros(p, ndb1);
        DDm(idn, :) = DD;
        AV(ip, :) = [zeros(p, nd-ndelta), DDm];
        Frtm1 = pinv(Frt);
        fstm = zeros(p);
        fstm(idn, idn) = Frtm1;
        FST(ip, :) = fstm;
        Kpm = zeros(p, nalpha);
        Kpm(idn, :) = K';
        KKP(ip, 1:nalpha) = Kpm;
    end
    ia = i * nalpha + 1:(i + 1) * nalpha;
    first = i * nkkp;
    iam = first + 1:first + nalpha;
    AAT(ia, :) = [zeros(nalpha, nd-ndelta), A(:, 1:ndb1)];
    PT(iam, 1:nalpha) = P;
    %
    % single collapse if possible
    %
    if collps == 0
        if rank(SQT(1:ndelta, 1:ndelta)) == ndelta
            LL = SQT(1:nd, 1:nd);
            hd = LL \ SQT(1:nd, nd+1:ndb1);
            Rm1 = LL \ eye(nd);
            M = Rm1 * Rm1';
            AHD = [M, hd(:, 1:nb1)];
            AA = A(:, 1:nd);
            AAT(ia, 1:nd) = AA;
            P = P + AA * M * AA';
            A = A(:, nd+1:ndb1) - AA * hd;
            AAT(ia, nd+1:ndb1) = A;
            PT(iam, 1:nalpha) = P;
            %
            % redefinition of SQT
            %
            if (nbeta > 0)
                nsqt = size(SQT, 1);
                SQT = SQT(nd+1:nsqt, nd+1:ndb1);
            else
                clear SQT
            end
            ndelta = 0;
            ndb1 = nb1;
            ndb = nbeta;
            iti = i;
            collps = 1;
        end
    end
    %   pause
end
%
% smoothing: first part (beta)
%
% computation of hat beta
%
ndb = nd + nbeta;
ndb1 = nd + nb1;
if nbeta > 0
    Rm1 = pinv(SQT(1:nbeta, 1:nbeta));
    hb = SQT(1:nbeta, 1:nbeta) \ SQT(1:nbeta, nbeta+1:nb1);
    Mb = Rm1 * Rm1';
else
    hb = [];
    Mb = [];
end
R = zeros(nalpha, nb1);
NN = zeros(nalpha, nalpha);
for i = n:-1:iti + 1
    %   i
    ip = (i - 1) * p + 1:i * p;
    ia = (i - 1) * nalpha + 1:i * nalpha;
    im = (i - 1) * mucd + 1:i * mucd;
    first = (i - 1) * nkkp;
    iam = first + 1:first + nalpha;
    iamm = first + 1:first + mucd;
    if mz > p
        ZZ = Z(ip, :);
    else
        ZZ = Z;
    end
    if mt > nalpha
        TT = T(ia, :);
    else
        TT = T;
    end
    if mg > p
        GG = G(ip, :);
    else
        GG = G;
    end
    if mh > nalpha
        HH = H(ia, :);
    else
        HH = H;
    end
    if mmc > mucd
        CC = C(im, :);
    else
        CC = C;
    end
    if md > mucd
        DD = D(im, :);
    else
        DD = D;
    end
    if mu > mucd
        UU = U(im, :);
    else
        UU = U;
    end
    Lm1 = FST(ip, :);
    V = AV(ip, nd+1:ndb1);
    K = KKP(ip, 1:nalpha)';
    A = AAT(ia, nd+1:ndb1);
    P = PT(iam, 1:nalpha);
    L = (TT - K * ZZ);
    Kjj = (CC * P * ZZ' + DD * GG') * Lm1';
    CKZ = (CC - Kjj * Lm1 * ZZ);
    DKG = (DD - Kjj * Lm1 * GG);
    Pjjm1 = CKZ * P * TT' + DKG * HH';
    Pjj = CKZ * P * CC' + DKG * DD';
    if isempty(UU)
        KKP(i*p, 1:mucd) = ((CC * A + Kjj * V + Pjjm1 * R) * [-hb', 1]')';
        DB = CC * A(:, 1:nbeta) + Kjj * V(:, 1:nbeta) + Pjjm1 * R(:, 1:nbeta);
    else
        KKP(i*p, 1:mucd) = ((CC * A + Kjj * V + Pjjm1 * R) * [-hb', 1]' + UU * hb)';
        DB = -UU + CC * A(:, 1:nbeta) + Kjj * V(:, 1:nbeta) + Pjjm1 * R(:, 1:nbeta);
    end
    PT(iamm, 1:mucd) = Pjj - Pjjm1 * NN * Pjjm1' + DB * Mb * DB';
    R = ZZ' * Lm1' * V + L' * R;
    NN = (Lm1 * ZZ)' * (Lm1 * ZZ) + L' * NN * L;
end
%
% smoothing: second part (delta,beta)
%
R = [zeros(nalpha, nd), R];
%
% computation of [hat delta' hat beta']'
%
if nd > 0
    ia = iti * nalpha + 1:(iti + 1) * nalpha;
    nhd = 1:nd;
    %
    % estimator of delta based on the initial stretch
    %
    Md = AHD(nhd, nhd);
    Sd1 = Md;
    hd = AHD(:, ndb1);
    if nbeta > 0
        cc1 = nd + 1:ndb;
        SA = Md * AAT(ia, nhd)';
        hd1 = hd - SA * R(:, ndb1);
        Ld = AHD(:, cc1) - SA * R(:, cc1);
        hd1 = hd1 - Ld * hb;
        hd = [hd1', hb']';
        Md1 = Md - SA * NN * SA';
        Md = [Md1 + Ld * Mb * Ld', -Ld * Mb; -Mb * Ld', Mb];
    else
        SA = Md * AAT(ia, nhd)';
        hd = hd - SA * R(:, ndb1);
        Md = Md - SA * NN * SA';
    end
else
    hd = hb;
    Md = Mb;
end
for i = iti:-1:1
    ip = (i - 1) * p + 1:i * p;
    ia = (i - 1) * nalpha + 1:i * nalpha;
    im = (i - 1) * mucd + 1:i * mucd;
    first = (i - 1) * nkkp;
    iam = first + 1:first + nalpha;
    iamm = first + 1:first + mucd;
    if mz > p
        ZZ = Z(ip, :);
    else
        ZZ = Z;
    end
    if mt > nalpha
        TT = T(ia, :);
    else
        TT = T;
    end
    if mg > p
        GG = G(ip, :);
    else
        GG = G;
    end
    if mh > nalpha
        HH = H(ia, :);
    else
        HH = H;
    end
    if mmc > mucd
        CC = C(im, :);
    else
        CC = C;
    end
    if md > mucd
        DD = D(im, :);
    else
        DD = D;
    end
    if mu > mucd
        UU = U(im, :);
    else
        UU = U;
    end
    Lm1 = FST(ip, :);
    V = AV(ip, :);
    K = KKP(ip, 1:nalpha)';
    A = AAT(ia, :);
    P = PT(iam, 1:nalpha);
    L = (TT - K * ZZ);
    Kjj = (CC * P * ZZ' + DD * GG') * Lm1';
    CKZ = (CC - Kjj * Lm1 * ZZ);
    DKG = (DD - Kjj * Lm1 * GG);
    Pjjm1 = CKZ * P * TT' + DKG * HH';
    Pjj = CKZ * P * CC' + DKG * DD';
    if isempty(UU)
        KKP(i*p, 1:mucd) = ((CC * A + Kjj * V + Pjjm1 * R) * [-hd', 1]')';
        DB = CC * A(:, 1:ndb) + Kjj * V(:, 1:ndb) + Pjjm1 * R(:, 1:ndb);
    else
        KKP(i*p, 1:mucd) = ((CC * A + Kjj * V + Pjjm1 * R) * [-hd', 1]' + UU * hb)';
        DB = [zeros(mmc, nd), -UU] + CC * A(:, 1:ndb) + Kjj * V(:, 1:ndb) + Pjjm1 * R(:, 1:ndb);
    end
    PT(iamm, 1:mucd) = Pjj - Pjjm1 * NN * Pjjm1' + DB * Md * DB';
    %
    % Adjustment in the covariance matrix due to the initial stretch.
    % The effect is due to delta only.
    %
    if nd > 0 && iti < n
        DBD = CC * A(:, 1:nd) + Kjj * V(:, 1:nd) + Pjjm1 * R(:, 1:nd);
        ia = i * nalpha + 1:(i + 1) * nalpha;
        AA = AAT(ia, :);
        CSBp = Pjjm1 * (R(:, 1:nd) + NN * AA(:, 1:nd)) * Sd1 * DBD';
        PT(iamm, 1:mucd) = PT(iamm, 1:mucd) - CSBp - CSBp';
    end
    R = ZZ' * Lm1' * V + L' * R;
    NN = (Lm1 * ZZ)' * (Lm1 * ZZ) + L' * NN * L;
end
for i = 1:n
    im = (i - 1) * mucd + 1:i * mucd;
    first = (i - 1) * nkkp;
    iamm = first + 1:first + mucd;
    KKP(i, 1:mucd) = KKP(i*p, 1:mucd);
    PT(im, 1:mucd) = PT(iamm, 1:mucd);
end
KKP = KKP(1:n, 1:mucd);
PT = PT(1:n*mucd, 1:mucd);