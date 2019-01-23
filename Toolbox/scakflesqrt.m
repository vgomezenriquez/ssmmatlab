function [e, f, hb, Mb, A, LP, qyy, R, olsres] = scakflesqrt(y, X, Z, G, W, T, H, ins, i, chb, icollps)
%
%
%       This function applies the square root version of the two stage
%       Kalman filter to the series y for prediction and likelihood
%       evaluation corresponding to the model
%
%       y_t = X_t*beta + Z_t*alpha_t + G_t*epsilon_t
%       alpha_{t+1}= W_t*beta + T_t*alpha_t + H_t*epsilon_t,
%
%       where epsilon_t is (0,sigma^2I),
%
%       with initial state
%
%       alpha_1= c + W_0*beta + a_1 + A_1*delta
%
%       where c is (0,Omega) and delta is (0,kI) (diffuse). A single
%       collapse is applied to get rid of the diffuse component. First, the
%       tskf-sribf if applied until collapsing. Then, the GLS model is
%       transformed into an OLS model by means of the tskf. The (E,e)
%       matrices of the OLS model are stored and, at the end of filtering,
%       the QR algorithm is applied to obtain the residuals and, if chb=1,
%       the regression estimates and their mse as well.
%
%       Input parameters:
%       y:     an (n x p) matrix of observations;
%       X    : an (n*p x nbeta) matrix containing the X_t matrices;
%              a  (p x nbeta) if it is time invariant;
%              it can be []
%       Z    : an (n*p x nalpha) matrix containing the Z_t matrices;
%              a  (p x nalpha) matrix if it is time invariant
%       G    : an (n*p x nepsilon) matrix containing the G_t matrices;
%              a  (p x nepsilon) matrix if it is time invariant
%       W    : an (n*nalpha x nbeta) matrix containing the W_t matrices;
%              an (nalpha x nbeta) matrix if it is time invariant;
%              it can be []
%       T    : an (n*nalpha x nalpha) matrix containing the T_t matrices;
%              an (nalpha x nalpha) matrix if it time invariant
%       H    : an (n*nalpha x nepsilon) matrix containing the H_t matrices;
%              an (nalpha x nepsilon) if it is time invariant
%       ins: an nalpha x (cc+cw0+ca1+cca1) matrix containing the initial
%            state information, according to array i below
%       i    : a  1 x 4 array containing 4 integers, i=[cc cw0 ca1 cca1],
%            where
%            cc   = nalpha if c is not missing (0 if c missing)
%            cw0  = number of columns in W_0 (0 if W_0 missing)
%            ca1  = 1 if a_1 is not missing (0 if a_1 missing)
%            cca1 = number of columns in A_1 (0 if A_1 missing)
%       chb= 1 compute hb and Mb
%            0 do not compute hb and Mb
%   icollps= an integer, corresponding to the observation number in which a
%            collapse takes place
%
%       Output parameters:
%       e   : residual vector (Q'_2*y)
%       f   : factor by which the residuals are to be multiplied
%             for minimization of the nonlinear sum of squares
%       hb  : the beta estimator
%       Mb  : the Mse of the beta estimator
%        A  : the estimated augmented state vector at the end of filtering
%       LP  : the square root of the Mse of A at the end of filtering
%      qyy  : Q'_1*y in the QR decomposition to estimate beta (chb=1)
%        R  : R in the QR decomposition to estimate beta      (chb=1)
%   olsres  : vector of OLS residuals                         (chb=1)
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
if nargin < 11
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
    LP = P;
    [nlp, ~] = size(LP);
else
    P = ins(:, 1:cc);
    %we use svd instead of Cholesky to obtain P^{1/2} because
    %P can be of less than full rank. If P is of full
    %rank, the orthogonal matrices U and V coincide.
    [U, SV, V] = svd(P);
    SSV = diag(sqrt(diag(SV)));
    LP = U * SSV * V'; % P^{1/2}
    [nlp, ~] = size(LP);
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
missingt = sum(sum(isnan(y))); %number of missing data
idx = 0; %number of rows nonzero in SQTi,SQT
collps = 0;
if nd == 0
    collps = 1;
end
SQT = zeros(ndb+p, ndb1);
SQTd = zeros(np-missingt-nd, nbeta+1);
f = 1;
fc = 0;
nmiss = 0;
missing = 0;
 
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
        %
        %tskf
        %
        V = [zeros(p-miss, ndelta), full(XX), YY] - ZZ * A;
        MM = [LP' * ZZ', LP' * TT'; GG', HH'];
        [~, RR] = qr(MM);
        [nzz, ~] = size(ZZ);
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
            Frt = RR(1:nzz, 1:nzz)';
            whK = RR(1:nzz, nzz+1:nzz+nlp)';
            f = f * prod(abs(diag(Frt)));
            [f, fc] = updatef(f, fc);
            LP = RR(nzz+1:nzz+nlp, nzz+1:nzz+nlp)';
            DD = Frt \ V;
            A = [zeros(nalpha, ndelta), -full(WW), zeros(nalpha, 1)] + TT * A + whK * DD;
        elseif singular == 1
            DD = zeros(p-miss, ndb1);
            MM = [LP' * TT'; HH'];
            [~, RR] = qr(MM);
            LP = RR(1:nlp, 1:nlp)';
            A = [zeros(nalpha, ndelta), -full(WW), zeros(nalpha, 1)] + TT * A;
            miss = p;
        elseif singular == 2
            error('innovations covariance singular in scakflesqrt')
        end
    elseif miss == p
        DD = zeros(p, ndb1);
        MM = [LP' * TT'; HH'];
        [~, RR] = qr(MM);
        LP = RR(1:nlp, 1:nlp)';
        A = [zeros(nalpha, ndelta), -full(WW), zeros(nalpha, 1)] + TT * A;
    else
        V = [zeros(p, ndelta), full(XX), YY] - ZZ * A;
        MM = [LP' * ZZ', LP' * TT'; GG', HH'];
        [~, RR] = qr(MM);
        [nzz, ~] = size(ZZ);
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
            Frt = RR(1:nzz, 1:nzz)';
            whK = RR(1:nzz, nzz+1:nzz+nlp)';
            f = f * prod(abs(diag(Frt)));
            [f, fc] = updatef(f, fc);
            LP = RR(nzz+1:nzz+nlp, nzz+1:nzz+nlp)';
            DD = Frt \ V;
            A = [zeros(nalpha, ndelta), -full(WW), zeros(nalpha, 1)] + TT * A + whK * DD;
        elseif singular == 1
            DD = zeros(p, ndb1);
            MM = [LP' * TT'; HH'];
            [~, RR] = qr(MM);
            LP = RR(1:nlp, 1:nlp)';
            A = [zeros(nalpha, ndelta), -full(WW), zeros(nalpha, 1)] + TT * A;
            miss = p;
        elseif singular == 2
            error('innovations covariance singular in scakflesqrt')
        end
    end
    %
    % update number of nonmissing and missing values
    %
    nmiss = nmiss + p - miss;
    missing = missing + miss;
    %
    % srif until collapse (matrix SQT). Then, GLS --> OLS (matrix SQTd)
    % ilb is the number of rows of SQT different from zero.
    %
    if ndelta > 0
        lb = min(idx, ndb);
        [Q, SQTp] = qr([DD(:, 1:ndb); SQT(1:lb, 1:ndb)]);
        SQTi = [SQTp, Q' * [DD(:, ndb1); SQT(1:lb, ndb1)]];
        idx = size(SQTi, 1);
        SQT(1:idx, :) = SQTi;
    elseif miss < p
        ndd = size(DD, 1);
        if ndd > 0
            SQTd(idx+1:idx+ndd, :) = DD;
            idx = idx + ndd;
        end
    end
    %
    % single collapse if possible
    %
    if (collps == 0) && (idx >= nd)
        LL = SQT(1:nd, 1:nd);
        if (rank(LL) == nd) && (i > icollps)
            Rm1 = LL \ eye(nd);
            hd = Rm1 * SQT(1:nd, nd+1:ndb1);
            AA = A(:, 1:nd);
            AAM = AA / LL;
            MM = [LP'; AAM'];
            [~, RR] = qr(MM);
            LP = RR(1:nlp, 1:nlp)';
            A = A(:, nd+1:ndb1) - AA * hd;
            f = f * abs(prod(diag(LL)));
            [f, fc] = updatef(f, fc);
            %
            % end of srif; from now on, store (E,e) (OLS matrices) in SQTd.
            % idx is the number of rows of SQT different from zero.
            %
            nidx = idx - nd;
            if nidx > 0  
                SQTd(1:idx-nd, :) = SQT(nd+1:idx, nd+1:ndb1);
                idx = nidx;
            else
                idx = 0;
            end
            clear SQT
            ndelta = 0;
            ndb1 = nb1;
            collps = 1;
        end
    end
end
if collps == 0 && idx >= nd
    error('No collapsing is possible in scakflesqtr')
end
f = (f^(1 / (nmiss - nd))) * (2^(fc / (nmiss - nd)));
%
% computation of residuals, hat beta and its Mse
%
qyy = [];
R = [];
%
%computation of OLS residuals added on 24-3-2017
%
olsres = [];
if nbeta > 0
    ms = size(SQTd, 2);
    if ms < nbeta
        error('Singular model in scakflesqrt')
    end
    [Q, R] = qr(SQTd(1:idx, 1:nbeta));
    qy = Q' * SQTd(1:idx, nbeta+1);
    qyy = qy(1:nbeta);
    e = qy(nbeta+1:idx);
    if chb == 1
        Mb = pinv(R(1:nbeta, :));
        hb = Mb * qy(1:nbeta);
        Mb = Mb * Mb';
        %OLS residuals
        olsres = SQTd(1:idx, nbeta+1) - SQTd(1:idx, 1:nbeta) * hb;
    else
        hb = [];
        Mb = [];
    end
else
    e = SQTd;
    hb = [];
    Mb = [];
end
