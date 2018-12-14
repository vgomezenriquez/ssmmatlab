function [e, f, hd, Md, A, P, ne] = tskfsribf(y, X, Z, G, W, T, H, ins, i, chb)
%
%
%        This function applies the tskf-sribf to the series y for
%        prediction and likelihood evaluation corresponding to the model
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
%       chb= 1 compute hb and Mb
%            0 do not compute hb and Mb
%
%
%        Output parameters:
%        e  : vector containing the stack of the standardized residuals
%        f  : factor by which the residuals are to be multiplied
%              for minimization of the nonlinear sum of squares
%        hd : the beta estimate
%        Md : the Mse of hd
%        A  : the estimated augmented state vector at the end of filtering
%        P  : the Mse of A at the end of filtering
%        ne : vector containing the stack of the observation numbers
%             corresponding to the standardized residuals
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
    disp('too many regressors in tskfsribf')
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
    disp('the number of columns in H and G should agree')
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
collps = 0;
if nd == 0
    collps = 1;
end
SQT = zeros(ndb+p, ndb1);
if (nbeta > 0)
    SQTd = zeros(nbeta+p, nb1);
end
f = 1;
fc = 0;
nomiss = 0;
nmiss = 0;
hd = [];
Md = [];
nmissy = sum(sum(isnan(y)));
e = zeros(np-nmissy, 1);
ne = e;
idx = 0;
idxi = 0;
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
    nmiss = nmiss + miss;
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
            f = f * prod(diag(Frt));
            [f, fc] = updatef(f, fc);
            DD = Frt \ V;
            K = (TT * P * ZZ' + HG) / F;
        elseif cp == size(F, 1)
            if ndb == 0
                DD = [];
            else
                DD = zeros(p, ndb1);
            end
            K = zeros(nalpha, p-miss);
        else
            error('singular matrix different from zero in tskfsribf')
        end
    elseif miss == p
        if ndb == 0
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
            f = f * prod(diag(Frt));
            [f, fc] = updatef(f, fc);
            DD = Frt \ V;
            K = (TT * P * ZZ' + HG) / F;
        elseif cp == size(F, 1)
            if ndb == 0
                DD = [];
            else
                DD = zeros(p, ndb1);
            end
            K = zeros(nalpha, p);
        else
            error('singular matrix different from zero in tskfsribf')
        end
    end
    if miss == p
        A = [zeros(nalpha, ndelta), -full(WW), zeros(nalpha, 1)] + TT * A;
        P = TT * P * TT' + HH2;
    else
        A = [zeros(nalpha, ndelta), -full(WW), zeros(nalpha, 1)] + TT * A + K * V;
        P = TT * P * (TT - K * ZZ)' + HH2 - HG * K';
    end
    %
    % update number of nonmissing values
    %
    nomiss = nomiss + p - miss;
    %
    % sribf
    %
    if ndelta > 0
        lb = min(idxi, ndb);
        [Q, SQTp] = qr([DD(:, 1:ndb); SQT(1:lb, 1:ndb)]);
        SQTi = [SQTp, Q' * [DD(:, ndb1); SQT(1:lb, ndb1)]];
        idxi = size(SQTi, 1);
        SQT(1:idxi, :) = SQTi;
    elseif nbeta > 0
        if (miss ~= p)
            [Q, SQTp] = qr([SQTd(1:ndb, 1:ndb); DD(:, 1:ndb)]);
            idp = size(SQTp, 1);
            SQTd(1:idp, :) = [SQTp, Q' * [SQTd(1:ndb, ndb1); DD(:, ndb1)]];
            if (nomiss > nd + nbeta)
                ets = SQTd(ndb+1:idp, ndb1);
                nets = size(ets, 1);
                e(idx+1:idx+nets) = ets;
                sec = repmat(i, size(ets));
                ne(idx+1:idx+nets) = sec;
                idx = idx + nets;
            end
        end
    else
        nets = size(DD, 1);
        if (nets > 0)
            e(idx+1:idx+nets) = DD;
            sec = repmat(i, size(DD));
            ne(idx+1:idx+nets) = sec;
            idx = idx + nets;
        end
    end
    %
    % single collapse if possible
    %
    if collps == 0 && idxi >= nd
        LL = SQT(1:nd, 1:nd);
        if rank(LL) == ndelta
            Rm1 = LL \ eye(nd);
            hd = Rm1 * SQT(1:nd, nd+1:ndb1);
            M = Rm1 * Rm1';
            AA = A(:, 1:nd);
            P = P + AA * M * AA';
            A = A(:, nd+1:ndb1) - AA * hd;
            f = f * abs(prod(diag(LL)));
            [f, fc] = updatef(f, fc);
            if (nomiss > ndb)
                ets = SQT(ndb+1:idxi, end);
                nets = size(ets, 1);
                e(idx+1:idx+nets) = ets;
                sec = repmat(i, size(ets));
                ne(idx+1:idx+nets) = sec;
                idx = idx + nets;
            end
            %
            % pass from SQT to SQTd after collapsing
            %
            nidx = idxi - nd;
            if nidx > 0 && nomiss >= ndb && nbeta > 0
                SQTd(1:idxi-nd, :) = SQT(nd+1:idxi, nd+1:ndb1);
            end
            ndelta = 0;
            ndb1 = nb1;
            ndb = nbeta;
            collps = 1;
            clear SQT
        end
    end
end
f = (f^(1 / (nomiss - nd))) * (2^(fc / (nomiss - nd)));
e = e(1:idx);
ne = ne(1:idx);
hd = [];
Md = [];
if nbeta > 0 && (chb == 1)
    U = SQTd(1:nbeta, 1:nbeta);
    invu = pinv(U);
    hd = invu * SQTd(1:nbeta, nb1);
    Md = invu * invu';
end
