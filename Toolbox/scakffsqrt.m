function [KKP, PT, hd, Md, initf, recrs, recr, srecr] = scakffsqrt(y, X, Z, G, W, T, H, ins, i, icollps)
%
%
%        This function applies the square root version of the two stage
%        Kalman filter to the series y corresponding to the model to obtain
%        the filtered estimator of the state vector, its Mse, and recursive
%        residuals.
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
%   icollps= an integer, corresponding to the observation number in which a
%                 collapse takes place
%
%
%        Output parameters:
%        KKP : an (n x nalpha) matrix containing the estimated x_{t|t}
%        PT  : an (n*nalpha x nalpha) matrix containing the
%              Mse of x_{t|t}
%        hd  : the beta estimate
%        Md  : the Mse of hd
%       initf: flag to indicate when regression parameters are identified
%       recrs: standardized recursive residuals
%        recr: recursive residuals
%       srecr: covariance matrices of recursive residuals
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
if nargin < 10
    icollps = 0;
else
    icollps = icollps - 1;
end
%
% check for inconsistencies
%
np = n * p;
nnalpha = n * nalpha;
if nbeta > np
    disp('too many regressors in scakfle2')
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
    LP = P;
    nlp = size(LP, 1);
else
    P = ins(:, 1:cc);
    %we use svd instead of Cholesky to obtain P^{1/2} because
    %P can be of less than full rank. If P is of full
    %rank, the orthogonal matrices U and V coincide.
    [U, SV, V] = svd(P);
    SSV = diag(sqrt(diag(SV)));
    LP = U * SSV * V'; % P^{1/2}
    nlp = size(LP, 1);
end %
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
nd = ndelta;
ndb = ndelta + nbeta;
ndb1 = ndb + 1;
nb1 = nbeta + 1;
KKP = zeros(n, nalpha);
PT = zeros(n*nalpha, nalpha);
SQT = zeros(ndb+p, ndb1);
if (nbeta > 0)
    SQTd = zeros(nbeta+p, nb1);
end
collps = 0;
if nd == 0
    collps = 1;
end
best = 0;
initf = 0;
qt = 0;
st2 = 0;
ndbc = ndb;
nmiss = 0;
hd = [];
Md = [];
nomiss = 0;
nmissy = sum(sum(isnan(y)));
idxi = 0;
idx1 = 0;
idxpm1 = 0;
idxpm2 = 0;
idxpnm1 = 0;
idxpnm2 = 0;
if (nmissy > 0)
    if (p == 1)
        recrs = zeros(n, 1);
        recr = recrs;
        srecr = recrs;
    else
        recrs = zeros(np, 2);
        recr = zeros(np, 2);
        srecr = zeros(n*p*(p + 1)/2, 2);
    end
else
    recrs = zeros(n, p);
    recr = recrs;
    srecr = zeros(n*p, p);
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
        end
        if ~isempty(ZZ)
            ZZ = ZZ(idn, :);
        end
        V = [zeros(p-miss, ndelta), full(XX), YY] - ZZ * A;
        mgg = size(GG, 2);
        MM = [LP' * ZZ', LP' * TT', LP'; GG', HH', zeros(mgg, nlp)];
        [~, RR] = qr(MM);
        nzz = size(ZZ, 1);
        %make square root of Sigma_t equal to Cholesky factor to stabilize the
        %residuals
        singular = 0;
        cont = 0;
        for ii = 1:nzz
            if RR(ii, ii) < 0
                RR(ii, :) = -RR(ii, :);
            end
            if abs(RR(ii, ii)) < eps
                cont = cont + 1;
            end
        end
        if cont == nzz
            singular = 1;
        elseif cont > 0
            singular = 2;
        end
        if singular == 0
            Frt = RR(1:nzz, 1:nzz)'; %Sigma_t^{1/2}
            whK = RR(1:nzz, nzz+1:nzz+nlp)'; %\tilde K_t
            whKf = RR(1:nzz, nzz+nlp+1:nzz+2*nlp)'; %\tilde Kf_t
            LPf = RR(nzz+1:nzz+nlp, nzz+nlp+1:nzz+2*nlp)'; %P_{t|t}^[1/2}
            DD = Frt \ V; %\tilde E_t
        elseif singular == 1
            MM = [LP' * TT'; HH'];
            [~, RR] = qr(MM);
            DD = zeros(nzz, ndb1);
            V = zeros(size(V));
            nmiss = nmiss + p;
            miss = p;
        elseif singular == 2
            error('innovations covariance singular in scakffsqrt')
        end
    elseif miss == p
        MM = [LP' * TT'; HH'];
        [~, RR] = qr(MM);
        DD = zeros(p, ndb1);
        V = zeros(size(ZZ,1),size(A,2));
    else
        V = [zeros(p, ndelta), full(XX), YY] - ZZ * A;
        mgg = size(GG, 2);
        MM = [LP' * ZZ', LP' * TT', LP'; GG', HH', zeros(mgg, nlp)];
        [~, RR] = qr(MM);
        nzz = size(ZZ, 1);
        %make square root of Sigma_t equal to Cholesky factor to stabilize the
        %residuals
        singular = 0;
        cont = 0;
        for ii = 1:nzz
            if RR(ii, ii) < 0
                RR(ii, :) = -RR(ii, :);
            end
            if abs(RR(ii, ii)) < eps
                cont = cont + 1;
            end
        end
        if cont == nzz
            singular = 1;
        elseif cont > 0
            singular = 2;
        end
        if singular == 0
            Frt = RR(1:nzz, 1:nzz)'; %Sigma_t^{1/2}
            whK = RR(1:nzz, nzz+1:nzz+nlp)'; %\tilde K_t
            whKf = RR(1:nzz, nzz+nlp+1:nzz+2*nlp)'; %\tilde Kf_t
            LPf = RR(nzz+1:nzz+nlp, nzz+nlp+1:nzz+2*nlp)'; %P_{t|t}^[1/2}
            DD = Frt \ V; %\tilde E_t
        elseif singular == 1
            MM = [LP' * TT'; HH'];
            [~, RR] = qr(MM);
            DD = zeros(p, ndb1);
            V = zeros(size(V));
            nmiss = nmiss + p;
            miss = p;
        elseif singular == 2
            error('innovations covariance singular in scakffsqrt')
        end
    end
    %
    % update number of nonmissing
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
            U = SQTd(1:nbeta, 1:nbeta);
            invu = pinv(U);
            hd = invu * SQTd(1:nbeta, nb1);
            Md = invu * invu';
            [Q, SQTp] = qr([SQTd(1:ndb, 1:ndb); DD(:, 1:ndb)]);
            idp = size(SQTp, 1);
            SQTd(1:idp, :) = [SQTp, Q' * [SQTd(1:ndb, ndb1); DD(:, ndb1)]];
            if (nomiss > nd + nbeta)
                et = V * [-hd', 1]'; %recursive residuals
                set = Frt * Frt' + V(:, 1:ndb) * Md * V(:, 1:ndb)'; %and their mse
                ets = SQTd(ndb+1:idp, ndb1); %standr. recursive residuals
                if (nmissy > 0)
                    if (p == 1)
                        recrs(idx1+1) = ets;
                        recr(idx1+1) = et;
                        srecr(idx1+1) = set;
                        idx1 = idx1 + 1;
                    else
                        nets = size(ets, 1);
                        sec = repmat(i, size(ets));
                        recrs(idxpm1+1:idxpm1+nets, :) = [ets, sec];
                        recr(idxpm1+1:idxpm1+nets, :) = [et, sec];
                        idxpm1 = idxpm1 + nets;
                        fst = vech(set);
                        nfst = size(fst, 1);
                        sec = repmat(i, size(fst));
                        srecr(idxpm2+1:idxpm2+nfst, :) = [fst, sec];
                        idxpm2 = idxpm2 + nfst;
                    end
                else
                    recrs(idxpnm1+1, :) = ets';
                    recr(idxpnm1+1, :) = et';
                    idxpnm1 = idxpnm1 + 1;
                    srecr(idxpnm2+1:idxpnm2+p, :) = set;
                    idxpnm2 = idxpnm2 + p;
                end
            end
            if best == 0
                if rank(SQTd(ndelta+1:ndb, ndelta+1:ndb)) == ndb - ndelta
                    best = 1;
                    initf = i + 1;
                end
            end
        end
    else
        ets = DD;
        if ~isempty(DD)
            if (nmissy > 0)
                if (p == 1)
                    recrs(idx+1) = DD;
                    recr(idx+1) = V;
                    srecr(idx+1) = Frt * Frt';
                    idx1 = idx1 + 1;
                else
                    nets = size(DD, 1);
                    sec = repmat(i, size(DD));
                    recrs(idxpm1+1:idxpm1+nets, :) = [DD, sec]; 
                    recr(idxpm1+1:idxpm1+nets, :) = [V, sec];
                    idxpm1 = idxpm1 + nets;
                    fst = vech(Frt*Frt');
                    nfst = size(fst, 1);
                    sec = repmat(i, size(fst));
                    srecr(idxpm2+1:idxpm2+nfst, :) = [fst, sec];
                    idxpm2 = idxpm2 + nfst;
                end
            else
                recrs(idxpnm1+1, :) = DD';
                recr(idxpnm1+1, :) = V';
                idxpnm1 = idxpnm1 + 1;
                srecr(idxpnm2+1:idxpnm2+p, :) = Frt * Frt';
                idxpnm2 = idxpnm2 + p;
            end
        end
    end
    %
    % measurement update
    %
    ia = (i - 1) * nalpha + 1:i * nalpha;
    if (collps == 1) && (nomiss > nd + nbeta)
        if (miss < p)
            deno = i * p - nmiss - ndbc;
            qt = qt + ets' * ets;
            st2 = qt / deno;
            if (nbeta > 0)
                KKP(i, :) = ((A + whKf * DD) * [-hd', 1]')';
                DB = A(:, 1:ndb) + whKf * DD(:, 1:ndb);
                PT(ia, :) = (LPf * LPf' + DB * Md * DB') * st2;
            else
                KKP(i, :) = (A + whKf * DD)';
                PT(ia, :) = (LPf * LPf') * st2;
            end
        elseif (miss == p)
            if (nbeta > 0)
                KKP(i, :) = (A * [-hd', 1]')';
            else
                KKP(i, :) = A';
            end
            PT(ia, :) = P * st2;
        end
    else
        KKP(i, :) = NaN;
        PT(ia, :) = NaN;
    end
    %
    % time update for tskf
    %
    if miss == p
        A = [zeros(nalpha, ndelta), -full(WW), zeros(nalpha, 1)] + TT * A; %x_{t+1|t}
        LP = RR(1:nlp, 1:nlp)'; %P_{t+1}^{1/2}
    else
        A = [zeros(nalpha, ndelta), -full(WW), zeros(nalpha, 1)] + TT * A + whK * DD; %x_{t+1|t}
        LP = RR(nzz+1:nzz+nlp, nzz+1:nzz+nlp)';  %P_{t+1}^{1/2}
    end
    
    %
    % single collapse if possible
    %
    if collps == 0 && idxi >= nd
        LL = SQT(1:nd, 1:nd);
        if (rank(LL) == ndelta) && (i > icollps)
            Rm1 = LL \ eye(nd);
            hd = Rm1 * SQT(1:nd, nd+1:ndb1);
            AA = A(:, 1:nd);
            AAM = AA / LL;
            MM = [LP'; AAM'];
            [~, RR] = qr(MM);
            LP = RR(1:nlp, 1:nlp)';
            A = A(:, nd+1:ndb1) - AA * hd;
            %
            % pass from SQT to SQTd after collapsing
            %
            nidx = idxi - nd;
            if nidx > 0 && nbeta > 0 
                SQTd(1:idxi-nd, :) = SQT(nd+1:idxi, nd+1:ndb1);
            end
            %
            % recursive residuals, if any
            %
            if (nomiss > ndb)
                %        ets=SQT(ndb+1:end,end);        %standardized recursive residuals
                ets = SQT(ndb+1:idxi, end); %standardized recursive residuals
                nets = size(ets, 1);
                nc = size(Frt, 2);
                %
                % equation (4.142), p. 294, of Multivariate time series..., Gomez
                % (2016), to obtain the covariance matrix of recursive residuals
                %
                AS = Q' * [Frt \ eye(nc); zeros(idxi-nc, nc)];
                nas = size(AS, 1);
                [~, Ri] = qr(AS(nas-nc+1:end, :));
                Sigmasmoh = Ri(nc-nets+1:end, nc-nets+1:end);
                Sigmasoh = pinv(Sigmasmoh);
                et = Sigmasoh * ets; %recursive residuals
                set = Sigmasoh * Sigmasoh'; %and their mse
                qt = qt + ets' * ets; %sum of stdr. recursive residuals
                if (nmissy > 0)
                    if (p == 1)
                        recrs(idx1+1) = ets;
                        recr(idx1+1) = et;
                        srecr(idx1+1) = set;
                        idx1 = idx1 + 1;
                    else
                        sec = repmat(i, size(ets));
                        recrs(idxpm1+1:idxpm1+nets, :) = [ets, sec];
                        recr(idxpm1+1:idxpm1+nets, :) = [et, sec];
                        idxpm1 = idxpm1 + nets;
                        fst = vech(set);
                        nfst = size(fst, 1);
                        sec = repmat(i, size(fst));
                        srecr(idxpm2+1:idxpm2+nfst, :) = [fst, sec];
                        idxpm2 = idxpm2 + nfst;
                    end
                else
                    if nets < p
                        ats = [ets; zeros(p-nets, 1)];
                        at = [et; zeros(p-nets, 1)];
                        sat = [set, zeros(nets, p-nets); zeros(p-nets, p)];
                    else
                        ats = ets;
                        at = et;
                        sat = set;
                    end
                    recrs(idxpnm1+1, :) = ats';
                    recr(idxpnm1+1, :) = at';
                    idxpnm1 = idxpnm1 + 1;
                    srecr(idxpnm2+1:idxpnm2+p, :) = sat;
                    idxpnm2 = idxpnm2 + p;
                end
            end
            clear SQT AS Ri
            ndelta = 0;
            ndb1 = nb1;
            ndb = nbeta;
            collps = 1;
        end
    end
end
if nbeta > 0
    if (miss ~= p)
        U = SQTd(1:nbeta, 1:nbeta);
        invu = pinv(U);
        hd = invu * SQTd(1:nbeta, nb1);  
        Md = invu * invu';
    end
end

if (idx1 > 0)
    recrs = recrs(1:idx1);
    recr = recr(1:idx1);
    srecr = srecr(1:idx1);
elseif (idxpm1 > 0)
    recrs = recrs(1:idxpm1, :);
    recr = recr(1:idxpm1, :);
    srecr = srecr(1:idxpm2, :);
else
    recrs = recrs(1:idxpnm1, :);
    recr = recr(1:idxpnm1, :);
    srecr = srecr(1:idxpnm2, :);
end
