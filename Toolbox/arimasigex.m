function outa = arimasigex(out, Ycomp)
%
% function to perform the canonical decomposition of an ARIMA model
% previously identified with function arimaestos.
%
% phi(B)*phi_s(B^s)*(delta*delta_s*y_t -mu) =
% th(B)*th_s(B^s)*a_t
%
%     INPUTS:
%     out    : a structure, output of function arimaestos
%     Ycomp  : a cell array, containing the assignment of each regression
%              variable to a component. Possible values are 'trend',
%              'seas','tran' and 'irreg'.
%
%  OUTPUTS:
%      outa  : a structure containing model information for the input
%              with fields:
%       title: a string with the name of series
%        orig: original series
%     Ycomp  : same as input
%          gft: flag for graphs, = 0, no graphs, =1 graphs
%          lam: flag for logs, = 0 logs, = 1, no logs
%          npr: number of forecasts
%        sconp: standard deviation of the innovations
%         conp: innovation variance
%         orig: original series
%            Y: array containing the regression variables
%      bg_year: initial year
%       bg_per: initial period
%        datei: date structure
%          str: structure that contains the ARIMA model information. It is
%               the output of function suvarmapqPQ.
%       compcd: structure that is output of function candec
%         strc: structure that contains information about the canonical
%               decomposition model in state space form. It is the output
%               of function sucdm.
%      compmat: structure containing the state space model matrices and
%               initial conditions for the state space model in which
%               the trend-cycle has been decomposed into a smooth trend
%               and a cycle
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

outa = [];

title = out.title;
outa.title = title;

if nargin < 2
    Ycomp = [];
end
outa.Ycomp = Ycomp;

if isfield(out, 'gft')
    gft = out.gft;
else
    gft = 0;
end
outa.gft = gft;

s = out.freq; %seasonal frequency
lam = out.model.lam; %parameter for logs (=0, logs)
if isfield(out.model, 'npr')
    npr = out.model.npr; %number of forecasts
else
    npr = 0;
end
outa.lam = lam;
outa.npr = npr;
dr = out.model.d; %regular differencing
ds = out.model.ds; %seasonal differencing
p = out.model.p; %regular AR
ps = out.model.ps; %seasonal AR
q = out.model.q; %regular MA
qs = out.model.qs; %seasonal MA
phif(:, :, 1) = 1;
Thf(:, :, 1) = 1;
thf(:, :, 1) = 1;
Phif(:, :, 1) = 1;
for i = 1:p
    phif(:, :, 1+i) = out.model.phi(i);
end
for i = 1:q
    thf(:, :, 1+i) = out.model.th(i);
end
for i = 1:ps
    Phif(:, :, 1+i) = out.model.phis(i);
end
for i = 1:qs
    Thf(:, :, 1+i) = out.model.ths(i);
end
sconp = out.model.resinf.sconp; %innovations standard error
Sigmaf = out.model.resinf.conp; %mse of the innovations
y = out.orig; %original series
yor = y;
if (lam == 0)
    y = log(yor);
else
    y = yor;
end
Y = out.model.Y; %matrix for regression variables
bg_year = out.nziyip(2); %initial year
bg_per = out.nziyip(3); %initial period
datei = cal(bg_year, bg_per, s);

outa.sconp = sconp;
outa.conp = Sigmaf;
outa.orig = yor;
outa.Y = Y;
outa.bg_year = bg_year;
outa.bg_per = bg_per;
outa.datei = datei;


%create model structure, str
[str, ~] = suvarmapqPQ(phif, thf, Phif, Thf, Sigmaf, s);
outa.str = str;

% set up trend-cycle and seasonal polynomials for the canonical
% decomposition
[phir, phis, thr, ths, phirst] = arima2rspol(phif, Phif, thf, Thf, s, dr, ds);

% perform canonical decomposition
[compcd, ierrcandec] = candec(phir, phis, thr, ths, phirst, s, dr, ds, sconp);
if ierrcandec > 0
    disp('canonical decomposition failed in arimasigex')
    return
end
sigma2i = compcd.itvar;
irreg = 1;
if sigma2i < 0
    disp('irregular spectrum negative')
    %  return
    disp('model is changed. sigma2i is made equal to zero.')
    disp('this is a provisional solution to the negative ')
    disp('irregular spectrum problem')
    compcd.itvar = 0.;
    irreg = 0;
end
outa.compcd = compcd;

% put the canonical decomposition model y_t = Y_t*beta + p_t + s_t + r_t +
% i_t into state space form (Akaike state space form for each component)
% create structure
[X, Z, G, W, T, H, ins, ii, strc, ~] = sucdm(compcd, y, Y, str, npr);
outa.strc = strc;
outa.compmat.X = X;
outa.compmat.Z = Z;
outa.compmat.G = G;
outa.compmat.W = W;
outa.compmat.T = T;
outa.compmat.H = H;
outa.compmat.ins = ins;
outa.compmat.ii = ii;


% Smoothing
cdomp = ones(3, 1);
if ~isempty(compcd.ptden)
    ptden = compcd.ptden;
    np = length(ptden);
    trendcycle = ones(1, np);
else
    trendcycle = [];
    cdomp(1) = 0;
end
if ~isempty(compcd.stden)
    stden = compcd.stden;
    ns = length(stden);
    seas = ones(1, ns) * 2;
else
    seas = [];
    cdomp(2) = 0;
end
if ~isempty(compcd.rt)
    rt = compcd.rt;
    nr = length(rt);
    tran = ones(1, nr) * 3;
else
    tran = [];
    cdomp(3) = 0;
end
compall = [trendcycle, seas, tran];
idxc = logical(cdomp);
compallc = {trendcycle, seas, tran};
ccomp = compallc(idxc);
lcomp = length(ccomp);
istoc = zeros(1, lcomp);
nalpha = size(T, 1);
c = zeros(lcomp, nalpha);
csa = zeros(1, nalpha);
for i = 1:lcomp
    d = ccomp{i};
    istoc(i) = d(1);
    k = find(d(1) == compall);
    l = Z(k);
    c(i, k:k+length(d)-1) = l;
    if ~isempty(seas) && d(1) ~= 2
        csa = csa + c(i, :);
    end
end
if irreg == 1
    c = [c; zeros(1, nalpha)];
    lcomp = lcomp + 1;
end
nepsilon = size(G, 2);
D = zeros(lcomp, nepsilon);
%incorporate irregular if present
if irreg == 1
    D(lcomp, nepsilon) = sconp;
end
%seasonally adjusted series
if ~isempty(seas)
    c = [c; csa];
    D = [D; sum(D)];
    lcomp = lcomp + 1;
end
c = logical(c);
mucd = lcomp;
[nx, mx] = size(X);
U = [];
% mucd,c,D,Z,pause
% stochastic components (no deterministic effects other than those present
% in the model definitions)
% Function smoothgen smooths a general vector:
% Y_t = U_t*beta + C_t*alpha_t + D_t*epsilon_t
% In this case, it is desired to smooth:
% Y_t = U_t*beta + C_t*alpha_t
% Hence, U_t = [], C_t=c and D_t=0
% For the irregular, U_t = [], C_t = 0 and D_t = (0,...,0,sigmaa)
% mucd: an integer, the dimension of Y_t
C = c;
[KKP, PT] = smoothgen(y, X, Z, G, W, T, H, ins, ii, mucd, U, C, D);
StochCc = KKP;
outa.StochCc = StochCc;
nKKP = size(KKP, 1);
StochSCc = zeros(size(KKP));
for i = 1:nKKP
    CC = PT((i - 1)*mucd+1:i*mucd, :);
    StochSCc(i, :) = sqrt(diag(CC)) * sconp;
end
outa.StochSCc = StochSCc;
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
    outa.oStochCc = oStochCc;
    outa.oStochSCc = oStochSCc;
end

%incorporate regression effects to the components
%
if ~isempty(Ycomp) %Y contains the regression variables
    mYcomp = size(Ycomp, 2);
    if mYcomp ~= mx
        disp('the number of elements in Ycomp should be equal to ')
        disp('the number of columns in Y')
        return
    end
    Ycompn = zeros(1, mYcomp);
    for i = 1:mYcomp
        if strcmp(Ycomp{i}, 'trend')
            Ycompn(i) = 1;
        elseif strcmp(Ycomp{i}, 'seas')
            Ycompn(i) = 2;
        elseif strcmp(Ycomp{i}, 'tran')
            Ycompn(i) = 3;
        elseif strcmp(Ycomp{i}, 'irreg')
            Ycompn(i) = 4;
        elseif isempty(Ycomp{i})
            Ycompn(i) = 0;
        else
            disp('the names in Ycomp are incorrect')
            return
        end
    end
    U = zeros(mucd*nx, mx);
    for i = 1:nx
        ip = (i - 1) * mucd;
        for j = 1:mx
            ic = Ycompn(j);
            if ismember(ic, istoc) || (ic == 4 && irreg == 1)
                contc = 0;
                for k = 1:4
                    if ismember(k, istoc) || (k == 4 && irreg == 1)
                        contc = contc + 1;
                    end
                    if k == ic
                        U(ip+contc, j) = X(i, j);
                        %       i,j,k, aa=U(ip+contc,j),pause
                    end
                end
            end
        end
        if ~isempty(seas)
            sac = zeros(1, mx);
            contc = 0;
            for k = 1:4
                if (ismember(k, istoc) || (k == 4 && irreg == 1)) && k ~= 2
                    contc = contc + 1;
                    sac = sac + U(ip+contc, :);
                elseif k == 2
                    contc = contc + 1;
                end
            end
            U(ip+mucd, :) = sac;
        end
    end
end

% Function smoothgen smooths a general vector:
% Y_t = U_t*beta + C_t*alpha_t + D_t*epsilon_t
% In this case, it is desired to smooth:
% Y_t = U_t*beta + C_t*alpha_t
% Hence, U_t = parts of X, C_t=c and D_t=0
% For the irregular, U_t = parts of X, C_t = 0 and D_t = (0,...,0,sigmaa)
% mucd: an integer, the dimension of Y_t
C = c;
[KKP, PT] = smoothgen(y, X, Z, G, W, T, H, ins, ii, mucd, U, C, D);
Cc = KKP;
outa.Cc = Cc;
nKKP = size(KKP, 1);
SCc = zeros(size(KKP));
for i = 1:nKKP
    CC = PT((i - 1)*mucd+1:i*mucd, :);
    SCc(i, :) = sqrt(diag(CC)) * sconp;
end
outa.SCc = SCc;
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
    outa.oCc = oCc;
    outa.oSCc = oSCc;
end

if gft == 1
    dbname = strrep(title, '_', '\_');
    names = char([dbname, ': Original series']);
    tsplot(yor, datei, names);
    pause
    close all
    if (lam == 0)
        names = char('Original series in logs');
        tsplot(y, datei, names);
        pause
        close all
    end
    %plot stochastic components
    compn = {'stoch. trend', 'stoch. seasonal', 'stoch. tran', 'stoch. irreg', ...
        'stoch. sa'};
    cont = 0;
    for i = 1:5
        ixc = ismember(i, istoc);
        if ixc > 0 && i == 1 && ismember(1, istoc)
            cont = cont + 1;
            names = char('original series', compn{i});
            tsplot([y, StochCc(:, 1)], datei, names);
            pause
            close all
            if lam == 0
                names = char('original series in levels ', [compn{i}, ...
                    ' in the original scale']);
                tsplot([yor, oStochCc(:, 1)], datei, names);
                pause
                close all
            end
        elseif ixc > 0 || (i == 4 && irreg == 1) || (i == 5 && ~isempty(seas))
            cont = cont + 1;
            names = compn{i};
            tsplot(StochCc(:, cont), datei, names);
            pause
            close all
            if lam == 0
                names = [compn{i}, ' in the original scale'];
                tsplot(oStochCc(:, cont), datei, names);
                pause
                close all
            end
            if (i == 5 && ~isempty(seas))
                names = char('original series', compn{i});
                tsplot([y, StochCc(:, cont)], datei, names);
                pause
                close all
                if lam == 0
                    names = char('original series in levels ', [compn{i}, ...
                        ' in the original scale']);
                    tsplot([yor, oStochCc(:, cont)], datei, names);
                    pause
                    close all
                end
            end
        end
    end
    %plot components with deterministic effects
    if ~isempty(Ycomp)
        compn = {'trend with det. components', 'seasonal with det. components', ...
            'tran with det. components', 'irreg with det. components', ...
            'sa with det. components'};
        cont = 0;
        for i = 1:5
            ixc = ismember(i, istoc);
            if ixc > 0 || (i == 4 && irreg == 1) || (i == 5 && ~isempty(seas))
                cont = cont + 1;
            end
            if ismember(i, Ycompn) || (i == 5 && ~isempty(seas))
                if ixc > 0 && i == 1 && ismember(1, istoc)
                    names = char('original series', compn{i});
                    tsplot([y, Cc(:, 1)], datei, names);
                    pause
                    close all
                    if lam == 0
                        names = char('original series in levels ', [compn{i}, ...
                            ' in the original scale']);
                        tsplot([yor, oCc(:, 1)], datei, names);
                        pause
                        close all
                    end
                elseif ixc > 0 || (i == 4 && irreg == 1) || (i == 5 && ~isempty(seas))
                    names = compn{i};
                    tsplot(Cc(:, cont), datei, names);
                    pause
                    close all
                    if lam == 0
                        names = [compn{i}, ' in the original scale'];
                        tsplot(oCc(:, cont), datei, names);
                        pause
                        close all
                    end
                    if (i == 5 && ~isempty(seas))
                        names = char('original series', compn{i});
                        tsplot([y, Cc(:, cont)], datei, names);
                        pause
                        close all
                        if lam == 0
                            names = char('original series in levels ', [compn{i}, ...
                                ' in the original scale']);
                            tsplot([yor, oCc(:, cont)], datei, names);
                            pause
                            close all
                        end
                    end
                end
            end
        end
    end
end
