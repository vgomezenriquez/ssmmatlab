function outa = usmestni(dbname, ser)
%
% function to identify, estimate and forecast a structural model for one
% series. The structural model is of the form:
%
%  y_t = Y_t*beta + p_t + s_t + u_t + v_t + e_t,
%
% where Y_t is a vector of regression variables, p_t is the trend, s_t is
% the seassonal, u_t is the cyclical, v_t is the AR, and e_t is the
% irregular component.  This function allows for other values of freq
% besides 1, 4 and 12.
% Note that freq is a field of structures comp and datei. It is also the
% number of seasons. Therefore, it determines the seasonal component.
%
%
%     INPUTS:
%     dbname : name of the series
%     ser    : a structure, containing the instructions for this function
%
%  OUTPUTS:
%      outa  : a structure containing model information for the input
%              with fields:
%       title: a string with the name of series
%      nziyip: a 1 x 3 array with number of obs., initial year, initial per.
%        freq: number of seasons
%        orig: original series
%       model: structre with model information. It contains the following
%              fields
%           lam: flag for logarithmic transformation, = 0, take logs, = 1,
%                do not take logs
%             X: X matrix in the state space form
%             Z: Z matrix in the state space form
%             G: G matrix in the state space form
%             W: W matrix in the state space form
%             T: T matrix in the state space form
%             H: H matrix in the state space form
%           ins: ins matrix for the initial conditions
%             i: i array for the initial conditions
%        resinf: structure containing residual information
%         sconp: residual standard error
%       StochCc: matrix containing the stochastic components
%      StochSCc: matrix containing the mse of the stochastic components
%      oStochCc: matrix containing the stochastic components in the
%                original scale
%     oStochSCc: matrix containing the mse of the stochastic components in
%                the original scale
%            Cc: matrix containing the components including deterministic
%                effects
%           SCc: matrix containing the mse of Cc
%           oCc: matrix containing the Cc in the original scale
%          oSCc: matrix containing the mse of the oCc
%           npr: number of forecasts
%            Xp: matrix containing the forecasts of X
%            Wp: matrix containing the forecasts of W
%           pry: forecasts
%          spry: mse of the forecasts
%          alpr: matrix containing the forecasts of the state vector
%         malpr: three dimensional array containing each of the covariance
%                matrices of alpr
%         salpr: matrix containing the mse of alpr
%          opry: forecasts in the original scale
%         ospry: mse of the forecasts in the original scale
%         oalpr: matrix containing the alpr in the original scale
%        osalpr: matrix containing the mse of oalpr
%         ser: the input structure
%      result: a structure containing estimation results. It has
%              the following fields:
%          xvf: array containing the estimated parameters
%           xf: array containing the fixed parameters
%            e: array containing the residuals
%           Ss: residual sum of squares
%           Ff: the product F'*F
%      sigma2c: standard error of the parameter concentrated out of the
%               likelihood
%         Pevf: prediction error variance
%        SPevf: square root of Pevf
%           tv: t-values of the estimated parameters
%           se: standard errors of the estimated parameters
%           .F: vector of nonlinear functions whose sum of squares is
%               minimized at the end of estimation
%            h: vector of regression parameters
%            M: mse of h
%            A: Augmented state vector at the end of filtering
%            P: mse matrix of A at the end of filtering
%       olsres: OLS residuals
%          tvr: t-values of the regression parameters
%          ser: standard errors of the regression parameters
%       ferror: flag for errors
%
% Copyright (c) 21 July 2015 by Victor Gomez
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

% profile on       %this is to optimize the code. It gives some information
%                  %about the performance of all of the functions.

tic;
clear x

yor = ser.yor;


serar = ser;
serar.npr = 0;
serar.gft = 0;
serar.pr = 0;
serar.usm = 1;
outar = arimaestni(dbname, serar);

outa.title = outar.title;
outa.nziyip = outar.nziyip;
outa.freq = outar.freq;
outa.orig = outar.orig;
yorm = outar.yor;
y = outar.y;
lam = outar.model.lam;
outa.model.lam = lam;
ser.lam = lam;
bg_year = ser.bg_year;
bg_per = ser.bg_per;
freq = ser.freq;
nmiss = 0;
if isfield(outar, 'model')
    if isfield(outar.model, 'nmiss')
        nmiss = outar.model.nmiss;
        idxn = outar.model.idxn;
        outa.model.nmiss = nmiss;
        outa.model.idxn = idxn;
    end
end


%default values for the program

usmdefval

outa.ser = ser;

%start of program

%files for output
if (pr == 1)
    pathc = pwd; %current directory
    if exist([pathc, filesep, 'results'], 'dir') ~= 7
        mkdir(pathc, 'results'); %create directory results
    end
    sname = ['results', filesep, dbname, '.txt'];
    fid = fopen(sname, 'w');
else
    fid = 1;
end

%structure for printing
clear inft
inft = minft(fid, fh, wd, nd, scale);


datei = cal(bg_year, bg_per, freq); %initial date and frequency.
outa.ser.datei = datei;

%copy npr in mpr and make npr zero for estimation
if npr > 0
    mpr = npr;
    npr = 0;
else
    mpr = 0;
end

yf = y; %filled in series if there are missing values
[~, mY] = size(Y);
nreg = mY;
if (nmiss > 0)
    outa.model.yf = yf;
    y(idxn) = NaN(size(idxn)); %series with missing values
end

%create structure and put model into state space form
[str, ferror] = suusm(comp, y, Y, npr);
if ferror > 0
    return
end
if ~isempty(W)
    if mpr > 0
        nalpha = size(str.T, 1);
        mny = size(y, 1);
        str.W = W(1:nalpha*mny, :);
        Wor = W;
    else
        str.W = W;
    end
else
    Wor = [];
end

if nlestim == 1
    %estimate model
    %
    [result, str] = usmestim(y, str);
    
    outa.result = result;
    
    %estimated and fixed parameters
    xvf = result.xvf;
    xf = result.xf;
else
    xvf = str.xv;
    xf = str.xf;
    %
    % get residuals and estimator
    %
    chb = 1;
    pvar = str.pvar;
    pfix = str.pfix;
    [F, e, h, M, Pevf, A, P, olsres] = smfun(xvf, y, s, pfix, pvar, xf, chb, str);
    nr = length(pvar);
    Ss = e' * e;
    Ff = F' * F;
    ne = length(e); %residual sum of squares
    % disp('concentrated parameter:')
    sigma2c = Ss / (ne - nr); %estimated sigma square
    % compute prediction error variance (finite sample)
    % disp('prediction error variance (finite sample)')
    Pevf = Pevf * sigma2c;
    
    result.xvf = xvf;
    result.xf = xf;
    result.e = e;
    result.Ss = Ss;
    result.Ff = Ff;
    result.sigma2c = sigma2c;
    result.Pevf = Pevf;
    SPevf = sqrt(diag(Pevf));
    result.SPevf = SPevf;
    result.tv = [];
    result.se = [];
    result.Cv = [];
    if ~isempty(str.X) || ~isempty(str.W)
        %t-values of regression parameters
        M = M * sigma2c;
        ser = sqrt(diag(M));
        tvr = h ./ ser;
        result.tvr = tvr;
        result.ser = ser;
    else
        result.tvr = [];
        result.ser = [];
    end
    result.F = F;
    result.h = h;
    result.M = M;
    result.A = A;
    result.P = P;
    result.olsres = olsres;
end

%t-values of varma estimated parameters are in result.tv
%t-values of estimated regression parameters are in result.tvr
%Note that the standard errors are divided by the concentrated parameter
%(sqrt(result.sigma2c))

%create estimated model
[X, Z, G, W, T, H, ins, ii, ~] = pr2usm(xvf, xf, str);

outa.model.X = X;
outa.model.Z = Z;
outa.model.G = G;
outa.model.W = W;
outa.model.T = T;
outa.model.H = H;
outa.model.ins = ins;
outa.model.i = ii;


if (nmiss > 0)
    %  disp('******** Computation of interpolated values   ********');
    %
    % Computation with function smoothgen.m
    %
    % Function smoothgen smooths a general vector:
    % Y_t = U_t*beta + C_t*alpha_t + D_t*epsilon_t
    % In this case, it is desired to smooth:
    % Y_t = Z*alpha_t + G*epsilon_t
    % Hence, U_t = X, C_t=Z and D_t=G
    if ~isempty(X) && ~isempty(W)
        [nx, mx] = size(X);
        [nw, mw] = size(W);
        XX = [X, zeros(nx, mw)];
        WW = [zeros(nw, mx), W];
    else
        XX = X;
        WW = W;
    end
    U = XX;
    mucd = 1.;
    C = Z;
    D = G;
    [yint, syint] = smoothgen(y, XX, Z, G, WW, T, H, ins, ii, mucd, U, C, D);
    outa.model.yinterp = yint;
    
    %obtain interpolated and original values
    Interp = yint(idxn, 1);
    sInterp = Interp;
    for i = 1:nmiss
        sInterp(i) = syint(idxn(i), 1);
    end
    sInterp = sqrt(sInterp*result.sigma2c);
    str.nmiss = nmiss;
    str.idxn = idxn;
    str.Interp = Interp;
    str.sInterp = sInterp;
    outa.model.Interp = Interp;
    outa.model.sInterp = sInterp;
    %obtain interpolations in the original scale using the log-normal
    %distribution
    if lam == 0
        oInterp = Interp;
        osInterp = sInterp;
        for i = 1:nmiss
            oInterp(i) = exp(Interp(i)+(sInterp(i)^2)/double(2.));
            osInterp(i) = exp(double(2.)*Interp(i)+sInterp(i)^2) ...
                * (exp(sInterp(i)^2) - double(1.));
        end
        str.oInterp = oInterp;
        str.osInterp = osInterp;
        outa.model.oInterp = oInterp;
        outa.model.osInterp = osInterp;
        yorint = yor;
        yorint(idxn) = oInterp;
    end
end


%residual diagnostics
e = result.e;
F = result.F;
Ss = e' * e;
Ff = F' * F;
ne = length(e); %residual sum of squares
Pevf = result.Pevf; %prediction error variance
% disp('standard error (finite sample)')
SPevf = result.SPevf;
sconp = sqrt(result.sigma2c); %standard error of concentrated parameter
ny = length(y);
pvar = str.pvar;
nr = length(pvar);
X = str.X;
nbeta = size(X, 2) + size(W, 2);
ndrs = ne + nbeta;
lagl = min(36, max([floor(.2*ny), 3 * freq, 10]));
%ols residuals can also be used for inference
% olsres=result.olsres;
if ~isempty(result.h) && olsres == 1
    e = result.olsres;
end
infr = rescomp(e, lagl, nr, Ss, Pevf, SPevf, Ff, ndrs, nbeta);

%computation of Pierce Qs value when there is seasonality
if (freq > 1)
    s = freq;
    ne = length(e);
    r = infr.r;
    orders = 1:floor(lagl/s);
    if isfield(ser.comp, 'seas')
        nr = 1;
    else
        nr = 0;
    end
    np = length(orders);
    Qs = zeros(np, 1);
    pvalQs = zeros(np, 1);
    dfQs = zeros(np, 1);
    for i = 1:np
        Qs(i) = sum((r(s:s:i*s).^2)./(ne - s:-s:ne - i * s)') * ne * (ne + 2);
        dfQs(i) = max(0, i-nr); %degrees of freedom
        if dfQs(i) > 0
            pvalQs(i) = 1 - gammp(dfQs(i)*.5, Qs(i)*.5);
        else
            pvalQs(i) = 1;
        end
    end;
    infr.Qs = Qs;
    infr.dfQs = dfQs;
    infr.pvalQs = pvalQs;
    %Bell's seasonal dummies F-test
    nb = length(e);
    Ys = seasdm(size(outa.orig, 1), datei);
    Yb = Ys(end-size(e, 1)+1:end, :);
    nYb = size(Yb, 2);
    [bb, ~, er] = bmols(e, Yb);
    Ybb = Yb * bb;
    F = (Ybb' * Ybb / nYb) / (er' * er / (nb - nYb));
    pf = fdis_cdf(F, nYb, nb-nYb);
    %p-value for the F-test
    if (pf < .5)
        pH = 2 * pf;
    else
        pH = 2 * (1 - pf);
    end
    infr.pFtest = pH;
    infr.F = F;
    infr.dfn = nYb;
    infr.dfd = nb - nYb;
    %end of Bell's F-test
end


outa.model.resinf = infr;

%compute standard errors and t-values for the regression part
if nbeta > 0
    g = result.h;
    M = result.M;
    seb = result.ser;
    tb = result.tvr; %standard errors and t-values
    %   Mbeta=[g seb tb];
    outa.model.hb = g;
    outa.model.Mb = M;
    outa.model.Y = Y;
    outa.model.seb = seb;
    outa.model.tb = tb;
end

%smoothing of all components except the irregular
%
%we could use function scakfs, but we would not obtain the mse if there are
%fixed components.
% [Xt,Pt]=scakfs(y,X,Z,G,W,T,H,ins,ii);
%vector g contains the estimates of the vector (delta',beta')'. Thus, the
%vector of regression estimates, hat(beta), is at the end of g.
% %example with constant slope
% trend=Xt(:,1)+X*g(end);
% outa.model.Xt=Xt;
% outa.model.Pt=Pt;

outa.model.sconp = sconp;
% Computation with function smoothgen. We get the mse too in all cases.
%
% stochastic components
cdomp = ones(5, 1);
if str.trend == 1 % stochastic level
    level = 1;
elseif str.trend == 2 % Butterworth tangent
    nlev = 2;
    level = ones(1, nlev);
else
    level = []; % no constant level
    cdomp(1) = 0;
end
if str.slope == 1 % stochastic slope
    slope = 2;
else
    slope = []; % no constant slope
    cdomp(2) = 0;
end
if str.seas > 0
    if str.seas == 4 % Butterworth tangent seasonality
        nseas = freq;
    else
        nseas = freq - 1; % stochastic dummy or trigonometric seasonality
    end
    seas = repmat(3, 1, nseas);
else
    seas = []; % no seasonals of fixed dummy seasonality
    cdomp(3) = 0;
end
if str.cycle ~= 0
    cycle = [4, 4];
else
    cycle = []; % no cycle
    cdomp(4) = 0;
end
if str.arp ~= 0
    nar = str.arp;
    ar = repmat(5, 1, nar);
else
    ar = []; % no AR component
    cdomp(5) = 0;
end
if str.irreg == 1
    irreg = 1;
else
    irreg = 0;
end
idxc = logical(cdomp);
compallc = {level, slope, seas, cycle, ar};
ccomp = compallc(idxc);
compall = [level, slope, seas, cycle, ar];
lcomp = length(ccomp);
istoc = zeros(1, lcomp);
nalpha = size(T, 1);
c = zeros(lcomp, nalpha);
for i = 1:lcomp
    d = ccomp{i};
    istoc(i) = d(1);
    k = find(d(1) == compall);
    l = Z(k);
    if d(1) == 2 %slope component
        l = 1;
    end
    c(i, k:k+length(d)-1) = l;
end
if irreg == 1
    c = [c; zeros(1, nalpha)];
    lcomp = lcomp + 1;
end
%  c,Z,,pause
% fixed components, if they exist, are in the order seasonal, slope and
% level.
cdomf = ones(3, 1);
cont = 0;
if str.seas == -1 % fixed dummy seasonality
    cont = cont + 1;
    nseas = freq;
    fseas = repmat(3, 1, nseas-1);
else
    cdomf(1) = 0;
    fseas = [];
end
if str.slope == -1 % fixed slope
    cont = cont + 1;
    fslope = 2;
else
    cdomf(2) = 0;
    fslope = [];
end
if str.trend == -1 % fixed level
    cont = cont + 1;
    flevel = 1;
else
    cdomf(3) = 0;
    flevel = [];
end
if (cont > 0)
    idxf = logical(cdomf);
    compallf = {fseas, fslope, flevel};
    ccompf = compallf(idxf);
    lcompf = length(ccompf);
    ifix = zeros(1, lcompf);
    for i = 1:lcompf
        ifix(i) = ccompf{i}(1);
    end
else
    lcompf = 0;
    U = [];
    ifix = [];
end
% istoc,ifix,pause
[nx, mx] = size(X);
if ~isempty(ifix)
    if ismember(1, ifix) && ismember(2, ifix)
        %level and slope are fixed. They count as one fixed component
        lcompf = lcompf - 1;
        ccompf{lcompf} = [1, 2];
        mucd = lcomp + lcompf;
        ifix = zeros(1, lcompf);
        for i = 1:lcompf
            ifix(i) = ccompf{i}(1);
        end
        % ifix,ccompf,lcompf,pause
    elseif ismember(2, ifix)
        %slope is fixed. It should be added to the level, that is stochastic
        mucd = lcomp + lcompf - 1;
    else
        mucd = lcomp + lcompf;
    end
    indc = zeros(1, 5);
    indf = zeros(1, 5);
    contc = 0;
    contf = 0;
    for i = 1:5
        if ismember(i, ifix)
            contc = contc + 1;
            contf = contf + 1;
            indf(contf) = 1;
            if ismember(i, istoc)
                indc(contc) = 1;
            end
        else
            contc = contc + 1;
            contf = contf + 1;
            if ismember(i, istoc)
                indc(contc) = 1;
            end
        end
    end
    % indf,pause
    if indf(1) == 0 && indf(2) == 1
        indf(1:4) = indf(2:5);
        aa = (ifix == 2);
        ifix(aa) = 1;
    end
    %  ifix, indc,indf,pause
    %  lcomp,lcompf,pause
    %regression part, matrix U
    UU = zeros(nx, mx, lcompf);
    cont = 0;
    for i = 1:lcompf
        idf = ccompf{i};
        lc = length(idf);
        UU(:, cont+1:cont+lc, lcompf-i+1) = X(:, cont+1:cont+lc);
        cont = cont + length(idf);
    end
    %   aa=UU(1:10,:,1),bb=UU(1:10,:,2),cc=X(1:10,:),pause
    U = zeros(mucd*nx, mx);
    for i = 1:nx
        ip = (i - 1) * mucd;
        cont = 0;
        contc = 0;
        for j = 1:5
            ix = max(indc(j), indf(j));
            if ix > 0
                contc = contc + 1;
            end
            if indf(j) == 1
                cont = cont + 1;
                U(ip+contc, :) = UU(i, :, cont);
            end
        end
        %    i, aa=U(ip+1:ip+mucd,:),pause
    end
    %stochastic part, redefine matrix C if necessary
    if (mucd > lcomp)
        cc = zeros(mucd, nalpha);
        cont = 0;
        contc = 0;
        for i = 1:5
            ix = max(indc(i), indf(i));
            if ix > 0
                contc = contc + 1;
            end
            if indc(i) == 1
                cont = cont + 1;
                cc(contc, :) = c(cont, :);
            end
        end
        c = cc;
    end
else
    mucd = lcomp;
end
c = logical(c);
% mucd,c,pause
% ifix,istoc,pause
% stochastic components (no deterministic effects other than those present
% in the model definitions)
% Function smoothgen smooths a general vector:
% Y_t = U_t*beta + C_t*alpha_t + D_t*epsilon_t
% In this case, it is desired to smooth:
% Y_t = U_t*beta + C_t*alpha_t
% Hence, U_t = parts of X, C_t=c and D_t=0
% For the irregular, U_t = parts of X, C_t = 0 and D_t = (0,...,0,1)
% mucd: an integer, the dimension of Y_t
C = c;
nepsilon = size(G, 2);
D = zeros(mucd, nepsilon);
%incorporate irregular if present
if irreg == 1
    D(mucd, nepsilon) = 1.;
end
if ~isempty(X) && ~isempty(W)
    [nx, mx] = size(X);
    [nw, mw] = size(W);
    nu = size(U, 1);
    XX = [X, zeros(nx, mw)];
    WW = [zeros(nw, mx), W];
    if (nu > 0)
        UU = [U, zeros(nu, mw)];
    else
        UU = U;
    end
else
    XX = X;
    WW = W;
    UU = U;
end
[KKP, PT] = smoothgen(y, XX, Z, G, WW, T, H, ins, ii, mucd, UU, C, D);
StochCc = KKP;
outa.model.StochCc = StochCc;
nKKP = size(KKP, 1);
StochSCc = zeros(size(KKP));
for i = 1:nKKP
    CC = PT((i - 1)*mucd+1:i*mucd, :);
    StochSCc(i, :) = sqrt(diag(CC)) * sconp;
end
outa.model.StochSCc = StochSCc;
%obtain component estimators in the original scale using the log-normal
%distribution
if lam == 0
    oStochCc = KKP;
    oStochSCc = StochSCc;
    for i = 1:mucd
        oStochCc(:, i) = exp(KKP(:, i)+(StochSCc(:, i).^2)./double(2.));
        oStochSCc(:, i) = exp(double(2.).*KKP(:, i)+StochSCc(:, i).^2) .* ...
            (exp(StochSCc(:, i).^2) - double(1.));
    end
    outa.model.oStochCc = oStochCc;
    outa.model.oStochSCc = oStochSCc;
end

%incorporate regression effects to the components
%
if ~isempty(Ycomp)
    if isempty(U)
        U = zeros(mucd*nx, mx);
    end
    for i = 1:nx
        ip = (i - 1) * mucd;
        conf = mx - nreg;
        for j = 1:nreg
            ic = Ycomp(j);
            if ismember(ic, ifix) || ismember(ic, istoc) || ic == 6
                contc = 0;
                for k = 1:6
                    if ismember(k, ifix) || ismember(k, istoc) || k == 6
                        contc = contc + 1;
                    end
                    if k == ic
                        U(ip+contc, conf+j) = X(i, conf+j);
                        %       i,j,k, aa=U(ip+contc,conf+j),pause
                    end
                end
            end
        end
    end
end

% Function smoothgen smooths a general vector:
% Y_t = U_t*beta + C_t*alpha_t + D_t*epsilon_t
% In this case, it is desired to smooth:
% Y_t = U_t*beta + C_t*alpha_t
% Hence, U_t = parts of X, C_t=c and D_t=0
% For the irregular, U_t = parts of X, C_t = 0 and D_t = (0,...,0,1)
% mucd: an integer, the dimension of Y_t
C = c;
nepsilon = size(G, 2);
D = zeros(mucd, nepsilon);
%incorporate irregular if present
if irreg == 1
    D(mucd, nepsilon) = 1.;
end
if ~isempty(X) && ~isempty(W)
    [nx, mx] = size(X);
    [nw, mw] = size(W);
    nu = size(U, 1);
    XX = [X, zeros(nx, mw)];
    WW = [zeros(nw, mx), W];
    if (nu > 0)
        UU = [U, zeros(nu, mw)];
    else
        UU = U;
    end
else
    XX = X;
    UU = U;
    WW = W;
end
[KKP, PT] = smoothgen(y, XX, Z, G, WW, T, H, ins, ii, mucd, UU, C, D);
outa.model.Cc = KKP;
nKKP = size(KKP, 1);
SCc = zeros(size(KKP));
for i = 1:nKKP
    CC = PT((i - 1)*mucd+1:i*mucd, :);
    SCc(i, :) = sqrt(diag(CC)) * sconp;
end
outa.model.SCc = SCc;
%obtain component estimators in the original scale using the log-normal
%distribution
if lam == 0
    oCc = KKP;
    oSCc = SCc;
    for i = 1:mucd
        oCc(:, i) = exp(KKP(:, i)+(SCc(:, i).^2)./double(2.));
        oSCc(:, i) = exp(double(2.).*KKP(:, i)+SCc(:, i).^2) .* (exp(SCc(:, i).^2) ...
            -double(1.));
    end
    outa.model.oCc = oCc;
    outa.model.oSCc = oSCc;
end

outa.model.Ycomp = Ycomp;
outa.model.istoc = istoc;
outa.model.ifix = ifix;
outa.model.irreg = irreg;

%compute forecasts
if mpr > 0
    %hb, Mb, A and P are in structure result. Here, hb is the vector of
    %regression estimates and Mb is the matrix of mse. A is the
    %estimated state vector, x_{t|t-1}, obtained with the Kalman filter at the
    %end of the sample and P is the matrix of standard errors.
    hb = result.h;
    Mb = result.M;
    A = result.A;
    P = result.P;
    %obtain forecasts of the series and the state vector
    npr = mpr;
    outa.model.npr = npr;
    [fstr, ~] = suusm(comp, y, Y, npr);
    Xp = fstr.X;
    Wp = Wor;
    if ~isempty(Xp)
        Xp = Xp(end-npr+1:end, :);
        outa.model.Xp = Xp;
    end
    if ~isempty(Wp)
        Wp = Wp(end-npr*nalpha+1:end, :);
        outa.model.Wp = Wp;
    end
    if ~isempty(Xp) && ~isempty(Wp)
        [nx, mx] = size(Xp);
        [nw, mw] = size(Wp);
        XXp = [Xp, zeros(nx, mw)];
        WWp = [zeros(nw, mx), Wp];
    else
        XXp = Xp;
        WWp = Wp;
    end
    m = 1; %number of series
    [pry, mypr, alpr, malpr] = ssmpred(npr, m, A, P, XXp, Z, G, WWp, T, H, hb, Mb);
    spry = zeros(m, npr);
    nalpha = size(T, 1);
    salpr = zeros(nalpha, npr);
    for i = 1:npr
        spry(:, i) = sqrt(diag(mypr(:, :, i))) * sconp;
        salpr(:, i) = sqrt(diag(malpr(:, :, i))) * sconp;
    end
    outa.model.pry = pry;
    outa.model.spry = spry;
    outa.model.alpr = alpr;
    outa.model.malpr = malpr;
    outa.model.salpr = salpr;
    %obtain forecasts in the original scale using the log-normal
    %distribution
    opry = pry;
    ospry = spry;
    oalpr = alpr;
    osalpr = salpr;
    if lam == 0
        for i = 1:npr
            opry(i) = exp(pry(i)+(spry(i)^2)/double(2.));
            ospry(i) = exp(double(2.)*pry(i)+spry(i)^2) * (exp(spry(i)^2) - double(1.));
            oalpr(i) = exp(alpr(i)+(salpr(i).^2)./double(2.));
            osalpr(i) = exp(double(2.).*alpr(i)+salpr(i).^2) .* (exp(salpr(i)^2) - double(1.));
        end
        outa.model.opry = opry;
        outa.model.ospry = ospry;
        outa.model.oalpr = oalpr;
        outa.model.osalpr = osalpr;
    end
end


%print series title and series
if pr == 1
    if freq == 1
        dateib = cal(1900, 1, 6); %print data in six columns
        inftb = inft;
        inftb.fh = 0; %without years and header
        fnameb = [dbname, ' (starting year is ', num2str(bg_year), ')'];
    else
        dateib = datei;
        inftb = inft;
        fnameb = dbname;
    end
    if nmiss > 0
        lamn = 1;
        prtser(fid, fnameb, yor, y, ny, dateib, inftb, lamn);
        fnameb = [fnameb, ' with missing values filled'];
    end
    %print estimation results
    printusmer(fid, datei, fnameb, yorm, y, ny, lam, str, result, mx, nbeta);
    %print residual diagnostics
    printres(fid, infr);
    
    %print forecasts
    if (npr > 0)
        fprintf(fid, '\n');
        clear in
        in.fid = fid;
        if lam == 1
            in.cnames = char('Obs.      ', 'Forecasts ', 'Std. Error');
            in.fmt = char('%5.0f', '%12.4f', '%12.4f');
            mprint([(ny + 1:ny + npr)', pry', spry'], in);
        else
            in.cnames = char('Obs.         ', 'Forecasts    ', 'Std. Error   ', ...
                'For. (levels)', 'Std. Error   ');
            in.fmt = char('%5.0f', '%12.4f', '%12.4f', '%19.4f', '%19.4f');
            mprint([(ny + 1:ny + npr)', pry', spry', opry', ospry'], in);
        end
        fprintf(fid, '\n');
    end
    
    %print estimated components
    lamn = 1;
    compn = {'trend', 'slope', 'seasonal', 'cycle', 'ar'};
    cont = 0;
    for i = 1:5
        ixf = ismember(i, ifix);
        ixc = ismember(i, istoc);
        ix = max(ixf, ixc);
        if ix > 0 && i == 1
            cont = cont + 1;
            names = char('trend');
            prtser(fid, names, KKP(:, 1), y, ny, dateib, inftb, lamn);
            if lam == 0
                names = char('trend in the original scale');
                prtser(fid, names, oCc(:, 1), y, ny, dateib, inftb, lamn);
            end
        elseif ix > 0
            cont = cont + 1;
            names = compn{i};
            prtser(fid, names, KKP(:, cont), y, ny, dateib, inftb, lamn);
            if lam == 0
                names = [compn{i}, ' in the original scale'];
                prtser(fid, names, oCc(:, cont), y, ny, dateib, inftb, lamn);
            end
        end
    end
end

%graphs
if gft >= 1
    pathc = pwd; %current directory
    if (gft > 1)
        if exist([pathc, filesep, 'graphs'], 'dir') ~= 7
            mkdir(pathc, 'graphs'); %create directory graphs
        end
    end
    if (gft == 3) || (gft == 1)
        gflag = 1;
    else
        gflag = 0;
    end
    if gft == 2
        f = figure('visible', 'off');
    else
        f = figure('visible', 'on');
    end
    dbname = strrep(dbname, '_', '\_');
    names = char('Original series');
    tsplot(yor, datei, names);
    if (gft > 1)
        saveas(f, [pathc, filesep, 'graphs', filesep, 'Original'], 'pdf')
    end
    if (gflag == 1)
        disp('strike any key when ready')
        pause
        close all
    end
    if (lam == 0)
        if gft == 2
            f = figure('visible', 'off');
        else
            f = figure('visible', 'on');
        end
        names = char('Original series in logs');
        tsplot(y, datei, names);
        if (gft > 1)
            saveas(f, [pathc, filesep, 'graphs', filesep, 'Originalinlogs'], 'pdf')
        end
        if (gflag == 1)
            disp('strike any key when ready')
            pause
            close all
        end
    end
    if freq > 1
        dsm = 1;
    else
        dsm = 0;
    end
    for dr = 0:1
        for ds = 0:dsm
            if gft == 2
                f = figure('visible', 'off');
            else
                f = figure('visible', 'on');
            end
            sacspacdif(yf, dbname, dr, ds, freq, lag, cw);
            if (gft > 1)
                saveas(f, [pathc, filesep, 'graphs', filesep, ['Original', num2str(dr), ...
                    num2str(ds)]], 'pdf')
            end
            if (gflag == 1)
                disp('strike any key when ready')
                pause
                close all
            end
        end
    end
    if (nmiss > 0)
        if gft == 2
            f = figure('visible', 'off');
        else
            f = figure('visible', 'on');
        end
        names = char('Interpolated series');
        tsplot(yint, datei, names);
        if (gft > 1)
            saveas(f, [pathc, filesep, 'graphs', filesep, 'Interp'], 'pdf')
        end
        if (gflag == 1)
            disp('strike any key when ready')
            pause
            close all
        end
        if (lam == 0)
            if gft == 2
                f = figure('visible', 'off');
            else
                f = figure('visible', 'on');
            end
            names = char('Interpolated series in levels');
            tsplot(yorint, datei, names);
            if (gft > 1)
                saveas(f, [pathc, filesep, 'graphs', filesep, 'Interpo'], 'pdf')
            end
            if (gflag == 1)
                disp('strike any key when ready')
                pause
                close all
            end
        end
    end
    nh = length(result.h);
    nregt = size(X, 2) + size(W, 2);
    ndif = nh - nregt;
    if (nh > nregt)
        names = char('Model Regression effects');
        if ~isempty(X)
            if gft == 2
                f = figure('visible', 'off');
            else
                f = figure('visible', 'on');
            end
            Xg = X(:, 1:ndif) * g(1:ndif);
            tsplot(Xg, datei, names);
            if (gft > 1)
                saveas(f, [pathc, filesep, 'graphs', filesep, 'Modelregef'], 'pdf')
            end
            if (gflag == 1)
                disp('strike any key when ready')
                pause
                close all
            end
        end
    end
    if (nreg > 0)
        if gft == 2
            f = figure('visible', 'off');
        else
            f = figure('visible', 'on');
        end
        names = char('Regression effects');
        Yg = Y * g(ndif+1:ndif+nreg);
        tsplot(Yg, datei, names);
        if (gft > 1)
            saveas(f, [pathc, filesep, 'graphs', filesep, 'Regef'], 'pdf')
        end
        if (gflag == 1)
            disp('strike any key when ready')
            pause
            close all
        end
    end
    %plot residual diagnostics
    if ~isempty(result.h) && olsres == 1
        names = 'OLS residuals';
    else
        names = 'residuals';
    end
    plotres([], [], [], [], [], 1.96, names, gft, 0, [], 0, [], infr, 1, 1);
    if (mpr > 0)
        %plot forecasts
        outp.pry = pry;
        outp.spry = spry;
        if (lam == 0)
            outp.opry = opry;
            outp.ospry = ospry;
        else
            outp.opry = [];
            outp.ospry = [];
        end
        outp.y = y;
        outp.yor = yor;
        outp.ny = ny;
        outp.npr = npr;
        outp.cw = cw;
        outp.tname = dbname;
        outp.lam = lam;
        outp.s = freq;
        outp.gft = gft;
        pfctsusm(outp);
    end
    %plot stochastic components
    compn = {'Stoch. trend', 'Stoch. slope', 'Stoch. seasonal', 'Stoch. cycle', ...
        'Stoch. ar', 'Stoch. irreg'};
    cont = 0;
    for i = 1:6
        ixf = ismember(i, ifix);
        ixc = ismember(i, istoc);
        ix = max(ixf, ixc);
        if ix > 0 && i == 1
            if gft == 2
                f = figure('visible', 'off');
            else
                f = figure('visible', 'on');
            end
            cont = cont + 1;
            names = char('original series', compn{i});
            tsplot([y, StochCc(:, 1)], datei, names);
            pname = strrep(compn{i}, '. ', '');
            if (gft > 1)
                saveas(f, [pathc, filesep, 'graphs', filesep, pname], 'pdf')
            end
            if (gflag == 1)
                disp('strike any key when ready')
                pause
                close all
            end
            if lam == 0
                if gft == 2
                    f = figure('visible', 'off');
                else
                    f = figure('visible', 'on');
                end
                names = char('original series in levels ', [compn{i}, ...
                    ' in the original scale']);
                tsplot([yor, oStochCc(:, 1)], datei, names);
                pname = strrep(compn{i}, '. ', '');
                if (gft > 1)
                    saveas(f, [pathc, filesep, 'graphs', filesep, [pname, 'o']], 'pdf')
                end
                if (gflag == 1)
                    disp('strike any key when ready')
                    pause
                    close all
                end
            end
        elseif ix > 0 || (i == 6 && irreg == 1)
            if gft == 2
                f = figure('visible', 'off');
            else
                f = figure('visible', 'on');
            end
            cont = cont + 1;
            names = compn{i};
            tsplot(StochCc(:, cont), datei, names);
            pname = strrep(compn{i}, '. ', '');
            if (gft > 1)
                saveas(f, [pathc, filesep, 'graphs', filesep, pname], 'pdf')
            end
            if (gflag == 1)
                disp('strike any key when ready')
                pause
                close all
            end
            if lam == 0 && i ~= 6
                if gft == 2
                    f = figure('visible', 'off');
                else
                    f = figure('visible', 'on');
                end
                names = [compn{i}, ' in the original scale'];
                tsplot(oStochCc(:, cont), datei, names);
                pname = strrep(compn{i}, '. ', '');
                if (gft > 1)
                    saveas(f, [pathc, filesep, 'graphs', filesep, [pname, 'o']], 'pdf')
                end
                if (gflag == 1)
                    disp('strike any key when ready')
                    pause
                    close all
                end
            end
        end
    end
    %plot components with deterministic effects
    if ~isempty(Ycomp)
        compn = {'Trend with det. components', 'Slope with det. components', ...
            'Seasonal with det. components', 'Cycle with det. components', ...
            'Ar with det. components', 'Irreg with det. components'};
        cont = 0;
        for i = 1:6
            if ismember(i, Ycomp)
                ixf = ismember(i, ifix);
                ixc = ismember(i, istoc);
                ix = max(ixf, ixc);
                if ix > 0 && i == 1
                    if gft == 2
                        f = figure('visible', 'off');
                    else
                        f = figure('visible', 'on');
                    end
                    cont = cont + 1;
                    names = char('original series', compn{i});
                    tsplot([y, KKP(:, 1)], datei, names);
                    pname = strrep(compn{i}, '. ', '');
                    pname = strrep(pname, ' with ', '');
                    if (gft > 1)
                        saveas(f, [pathc, filesep, 'graphs', filesep, pname], 'pdf')
                    end
                    if (gflag == 1)
                        disp('strike any key when ready')
                        pause
                        close all
                    end
                    if lam == 0
                        if gft == 2
                            f = figure('visible', 'off');
                        else
                            f = figure('visible', 'on');
                        end
                        names = char('original series in levels ', [compn{i}, ...
                            ' in the original scale']);
                        tsplot([yor, oCc(:, 1)], datei, names);
                        pname = strrep(compn{i}, '. ', '');
                        pname = strrep(pname, ' with ', '');
                        if (gft > 1)
                            saveas(f, [pathc, filesep, 'graphs', filesep, [pname, 'o']], 'pdf')
                        end
                        if (gflag == 1)
                            disp('strike any key when ready')
                            pause
                            close all
                        end
                    end
                elseif ix > 0 || (i == 6 && irreg == 1)
                    if gft == 2
                        f = figure('visible', 'off');
                    else
                        f = figure('visible', 'on');
                    end
                    cont = cont + 1;
                    names = compn{i};
                    tsplot(KKP(:, cont), datei, names);
                    pname = strrep(compn{i}, '. ', '');
                    pname = strrep(pname, ' with ', '');
                    if (gft > 1)
                        saveas(f, [pathc, filesep, 'graphs', filesep, pname], 'pdf')
                    end
                    if (gflag == 1)
                        disp('strike any key when ready')
                        pause
                        close all
                    end
                    if lam == 0 && i ~= 6
                        if gft == 2
                            f = figure('visible', 'off');
                        else
                            f = figure('visible', 'on');
                        end
                        names = [compn{i}, ' in the original scale'];
                        tsplot(oCc(:, cont), datei, names);
                        pname = strrep(compn{i}, '. ', '');
                        pname = strrep(pname, ' with ', '');
                        if (gft > 1)
                            saveas(f, [pathc, filesep, 'graphs', filesep, [pname, 'o']], 'pdf')
                        end
                        if (gflag == 1)
                            disp('strike any key when ready')
                            pause
                            close all
                        end
                    end
                end
            end
        end
    end
end

% profile off
