function [KKP, PT, hd, Md] = scakfs(y, X, Z, G, W, T, H, ins, i)
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
%
%
%        Output parameters:
%        KKP : an (n x nalpha) matrix containing the estimated x_{t|n}
%        PT  : an (n*nalpha x nalpha) matrix containing the
%              Mse of x_{t|n}
%        hd  : the (delta',beta')' estimate
%        Md  : the Mse of hd
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
if nbeta > np
    disp('too many regressors in scakfs')
    return
end
if mi ~= nalpha
    disp('the number of rows in ins is incorrect')
    return
end
if mx > p || mz > p || mg > p || mw > nalpha ...
        || mt > nalpha || mh > nalpha
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
end
if nw ~= nx && nw ~= 0 && nx ~= 0
    disp('the number of columns in W and X should agree')
    return
end
if nt ~= nalpha
    disp('the number of columns in T and Z should agree')
    return
end
if nh ~= neps
    disp('the number of columns in G and H should agree')
    return
end
if mi ~= nalpha
    disp('the number of rows in ins and the number of columns in Z should agree')
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
nd = ndelta;
ndb = ndelta + nbeta;
ndb1 = ndb + 1;
nb1 = nbeta + 1;
AV = sparse(n*p, ndb1);
AHD = zeros(nd, ndb1);
FST = sparse(n*p, p);
KKP = zeros(n*p, nalpha);
AAT = sparse((n + 1)*nalpha, ndb1);
PT = zeros((n + 1)*nalpha, nalpha);
SQT = zeros(ndb+p, ndb1);
AAT(1:nalpha, :) = A;
PT(1:nalpha, :) = P;
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
        elseif cp == size(F, 1) %F is zero
            if ndelta == 0
                DD = [];
            else
                DD = zeros(p, ndb1);
            end
            K = zeros(nalpha, p-miss);
        else
            error('singular matrix different from zero in scakfle2')
        end
    elseif miss == p
        if ndelta == 0
            DD = [];
        else
            DD = zeros(p, ndb1);
        end
        K = zeros(nalpha, p);
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
        elseif cp == size(F, 1) %F is zero
            if ndelta == 0
                DD = [];
            else
                DD = zeros(p, ndb1);
            end
            K = zeros(nalpha, p);
        else
            error('singular matrix different from zero in scakfle2')
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
        KKP(ip, :) = zeros(p, nalpha);
    else
        A = [zeros(nalpha, ndelta), -full(WW), zeros(nalpha, 1)] + TT * A + K * V;
        P = TT * P * (TT - K * ZZ)' + HH2 - HG * K';
        if miss == 0
            idn = 1:p;
        end
        DDm = zeros(p, ndb1);
        DDm(idn, :) = DD;
        AV(ip, :) = [zeros(p, nd-ndelta), DDm];
        Frtm1 = Frt \ eye(p-miss);
        fstm = zeros(p);
        fstm(idn, idn) = Frtm1;
        FST(ip, :) = fstm;
        Kpm = zeros(p, nalpha);
        Kpm(idn, :) = K';
        KKP(ip, :) = Kpm;
    end
    ia = i * nalpha + 1:(i + 1) * nalpha;
    AAT(ia, :) = [zeros(nalpha, nd-ndelta), A(:, 1:ndb1)];
    PT(ia, :) = P;
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
            PT(ia, :) = P;
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
    hb = Rm1 * SQT(1:nbeta, nbeta+1:nb1);
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
    Lm1 = FST(ip, :);
    V = AV(ip, nd+1:ndb1);
    K = KKP(ip, :)';
    A = AAT(ia, nd+1:ndb1);
    P = PT(ia, :);
    L = (TT - K * ZZ);
    R = ZZ' * Lm1' * V + L' * R;
    NN = (Lm1 * ZZ)' * (Lm1 * ZZ) + L' * NN * L;
    KKP(i*p, :) = ((A + P * R) * [-hb', 1]')';
    DB = A(:, 1:nbeta) + P * R(:, 1:nbeta);
    PT(ia, :) = P - P * NN * P + DB * Mb * DB';
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
    %   i
    ip = (i - 1) * p + 1:i * p;
    ia = (i - 1) * nalpha + 1:i * nalpha;
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
    Lm1 = FST(ip, :);
    V = AV(ip, :);
    K = KKP(ip, :)';
    A = AAT(ia, :);
    P = PT(ia, :);
    L = (TT - K * ZZ);
    R = ZZ' * Lm1' * V + L' * R;
    NN = (Lm1 * ZZ)' * (Lm1 * ZZ) + L' * NN * L;
    KKP(i*p, :) = ((A + P * R) * [-hd', 1]')';
    DB = A(:, 1:ndb) + P * R(:, 1:ndb);
    PT(ia, :) = P - P * NN * P + DB * Md * DB';
    %
    % Adjustment in the covariance matrix due to the initial stretch.
    % The effect is due to delta only.
    %
    if nd > 0 && iti < n
        DBD = A(:, 1:nd) + P * R(:, 1:nd);
        CSBp = P * (R(:, 1:nd) + NN * A(:, 1:nd)) * Sd1 * DBD';
        PT(ia, :) = PT(ia, :) - CSBp - CSBp';
    end
end
PT = PT(1:n*nalpha, :);
for i = 1:n
    KKP(i, :) = KKP(i*p, :);
end
KKP = KKP(1:n, :);
