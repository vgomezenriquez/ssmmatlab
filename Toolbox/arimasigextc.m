function outb = arimasigextc(outa, comp, filter)
%
% function to perform the decomposition of the trend-cycle component of a
% canonical decomposition of an ARIMA model
% previously identified with function arimaestos.
%
% phi(B)*phi_s(B^s)*(delta*delta_s*y_t -mu) =
% th(B)*th_s(B^s)*a_t
%
% the decomposition of the trend-cycle, p_t, is of the form p_t = sp_t +
% c_t, where sp_t is a (smooth) trend and c_t is a (smooth) cycle.
%
%     INPUTS:
%     out    : a structure, output of function arimaestos
%     Ycomp  : a cell array, containing the assignment of each regression
%              variable to a component. Possible values are 'trend',
%              'seas','tran' and 'irreg'.
%     filter : flag for the filter to be applied to the trend-cycle
%              component of the canonical decomposition. Possible values
%              are 'lp' (low-pass) and 'bp' (band-pass).
%
%  OUTPUTS:
%      outb  : a structure containing model information for the input
%              with fields:
%           strc: structure that contains information about the model y_t =
%                 Y_t*beta + sp_t + c_t + s_t + r_t + i_t in Akaike state
%                 space form, where sp_t is the (smooth) trend and c_t is
%                 the (smooth) cycle,  both obtained from the previous p_t
%                 by application of the low pass or band pass filter. It is
%                 the output of function sucdmpbst (low pass) or sucdmpbp
%                 (band pass)
%      StochCctc: matrix containing the stochastic components
%     StochSCctc: matrix containing the mse of the stochastic components
%     oStochCctc: matrix containing the stochastic components in the
%                 original scale
%    oStochSCctc: matrix containing the mse of the stochastic components in
%                 the original scale
%           Cctc: matrix containing the components including deterministic
%                 effects
%          SCctc: matrix containing the mse of Cctc
%          oCctc: matrix containing the Cctc in the original scale
%         oSCctc: matrix containing the mse of the oCctc
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

outb = [];

title = outa.title;
sconp = outa.sconp;
lam = outa.lam;
% conp=outa.conp;
y = outa.orig;
yor = y;
if (lam == 0)
    y = log(yor);
else
    y = yor;
end
Y = outa.Y;
npr = outa.npr;
% bg_year=outa.bg_year;
% bg_per=outa.bg_per;
datei = outa.datei;

str = outa.str;
compcd = outa.compcd;
Ycomp = outa.Ycomp;
gft = outa.gft;
irreg = 1;
if compcd.itvar == 0
    irreg = 0;
end

% put the model y_t = Y_t*beta + sp_t + c_t + s_t + r_t + i_t into Akaike
% state space form, where sp_t is the (smooth) trend and c_t is the
% (smooth) cycle,  both obtained from the previous p_t by application of
% the low pass or band pass filter.
if filter == 'lp'
    [X, Z, G, W, T, H, ins, ii, strc, ~] = sucdmpbst(compcd, comp, y, Y, str, npr);
elseif filter == 'bp'
    [X, Z, G, W, T, H, ins, ii, strc, ~] = sucdmpbp(compcd, comp, y, Y, str, npr);
else
    disp('filter must be equal to ''lp'' (low-pass) or ''bp'' (band-pass)')
    disp('in arimasigextc')
    return
end

outb.strc = strc;


% Smoothing using function smoothgen
compcd = strc.compcd;
compf = strc.compf;
cdomp = ones(4, 1);
if ~isempty(compcd.stden)
    stden = compcd.stden;
    ns = length(stden);
    seas = ones(1, ns) * 3;
else
    seas = [];
    ns = 0;
    cdomp(3) = 0;
end
if ~isempty(compcd.rt)
    rt = compcd.rt;
    nr = length(rt);
    tran = ones(1, nr) * 4;
else
    tran = [];
    nr = 0;
    cdomp(4) = 0;
end
if ~isempty(compcd.ptden) && ~isempty(compf.den)
    ptden = compcd.ptden;
    den = compf.den;
    nsp = length(den) + length(ptden) - 1;
    trend = ones(1, nsp) * 1;
    nc = length(Z) - (nsp + ns + nr);
    cycle = ones(1, nc) * 2;
else
    trend = [];
    cycle = [];
    cdomp(1) = 0;
    cdomp(2) = 0;
end
compall = [trend, cycle, seas, tran];
idxc = logical(cdomp);
compallc = {trend, cycle, seas, tran};
ccomp = compallc(idxc);
lcomp = length(ccomp);
istoc = zeros(1, lcomp);
nalpha = size(T, 1);
c = zeros(lcomp, nalpha);
for i = 1:lcomp
    d = ccomp{i};
    istoc(i) = d(1);
    k = find(d(1) == compall);
    l = Z(k);
    c(i, k:k+length(d)-1) = l;
end
if irreg == 1
    c = [c; zeros(1, nalpha)];
    lcomp = lcomp + 1;
end
c = logical(c);
mucd = lcomp;
[nx, mx] = size(X);
U = [];
% mucd,c,pause
% stochastic components (no deterministic effects other than those present
% in the model definitions)
% Function smoothgen smooths a general vector:
% Y_t = U_t*beta + C_t*alpha_t + D_t*epsilon_t
% In this case, it is desired to smooth:
% Y_t = U_t*beta + C_t*alpha_t
% Hence, U_t = parts of X, C_t=c and D_t=0
% For the irregular, U_t = parts of X, C_t = 0 and D_t = (0,...,0,sigmaa)
% mucd: an integer, the dimension of Y_t
C = c;
nepsilon = size(G, 2);
D = zeros(mucd, nepsilon);
%incorporate irregular if present
if irreg == 1
    D(mucd, nepsilon) = sconp;
end
[KKP, PT] = smoothgen(y, X, Z, G, W, T, H, ins, ii, mucd, U, C, D);
StochCctc = KKP;
outb.StochCctc = StochCctc;
nKKP = size(KKP, 1);
StochSCctc = zeros(size(KKP));
for i = 1:nKKP
    CC = PT((i - 1)*mucd+1:i*mucd, :);
    StochSCctc(i, :) = sqrt(diag(CC)) * sconp;
end
outb.StochSCctc = StochSCctc;

%obtain component estimators in the original scale using the log-normal
%distribution
if lam == 0
    oStochCctc = KKP;
    oStochSCctc = StochSCctc;
    for i = 1:mucd
        oStochCctc(:, i) = exp(KKP(:, i)+(StochSCctc(:, i).^2)./double(2.));
        oStochSCctc(:, i) = exp(double(2.).*KKP(:, i)+StochSCctc(:, i).^2) .* ...
            (exp(StochSCctc(:, i).^2) - double(1.));
    end
    outb.oStochCctc = oStochCctc;
    outb.oStochSCctc = oStochSCctc;
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
        elseif strcmp(Ycomp{i}, 'cycle')
            Ycompn(i) = 2;
        elseif strcmp(Ycomp{i}, 'seas')
            Ycompn(i) = 3;
        elseif strcmp(Ycomp{i}, 'tran')
            Ycompn(i) = 4;
        elseif strcmp(Ycomp{i}, 'irreg')
            Ycompn(i) = 5;
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
            if ismember(ic, istoc) || (ic == 5 && irreg == 1)
                contc = 0;
                for k = 1:4
                    if ismember(k, istoc) || (k == 5 && irreg == 1)
                        contc = contc + 1;
                    end
                    if k == ic
                        U(ip+contc, j) = X(i, j);
                        %       i,j,k, aa=U(ip+contc,j),pause
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
% For the irregular, U_t = parts of X, C_t = 0 and D_t = (0,...,0,sconp)
% mucd: an integer, the dimension of Y_t
C = c;
[KKP, PT] = smoothgen(y, X, Z, G, W, T, H, ins, ii, mucd, U, C, D);
Cctc = KKP;
outb.Cctc = Cctc;
nKKP = size(KKP, 1);
SCctc = zeros(size(KKP));
for i = 1:nKKP
    CC = PT((i - 1)*mucd+1:i*mucd, :);
    SCctc(i, :) = sqrt(diag(CC)) * sconp;
end
outb.SCctc = SCctc;
%obtain component estimators in the original scale using the log-normal
%distribution
if lam == 0
    oCctc = KKP;
    oSCctc = SCctc;
    for i = 1:mucd
        oCctc(:, i) = exp(KKP(:, i)+(SCctc(:, i).^2)./double(2.));
        oSCctc(:, i) = exp(double(2.).*KKP(:, i)+SCctc(:, i).^2) .* ...
            (exp(SCctc(:, i).^2) - double(1.));
    end
    outa.oCctc = oCctc;
    outa.oSCctc = oSCctc;
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
    compn = {'stoch. trend', 'stoch. cycle', 'stoch. seasonal', 'stoch. tran', ...
        'stoch. irreg'};
    cont = 0;
    for i = 1:5
        ixc = ismember(i, istoc);
        if ixc > 0 && i == 1 && ismember(1, istoc)
            cont = cont + 1;
            names = char('original series', compn{i});
            tsplot([y, StochCctc(:, 1)], datei, names);
            pause
            close all
            if lam == 0
                names = char('original series in levels ', [compn{i}, ...
                    ' in the original scale']);
                tsplot([yor, oStochCctc(:, 1)], datei, names);
                pause
                close all
            end
        elseif ixc > 0 || (i == 5 && irreg == 1)
            cont = cont + 1;
            names = compn{i};
            tsplot(StochCctc(:, cont), datei, names);
            pause
            close all
            if lam == 0
                names = [compn{i}, ' in the original scale'];
                tsplot(oStochCctc(:, cont), datei, names);
                pause
                close all
            end
        end
    end
    %plot components with deterministic effects
    if ~isempty(Ycomp)
        compn = {'trend with det. components', 'cycle with det. components', ...
            'seasonal with det. components', ...
            'tran with det. components', 'irreg with det. components'};
        cont = 0;
        for i = 1:5
            ixc = ismember(i, istoc);
            if ixc > 0 || (i == 5 && irreg == 1)
                cont = cont + 1;
            end
            if ismember(i, Ycompn)
                if ixc > 0 && i == 1 && ismember(1, istoc)
                    names = char('original series', compn{i});
                    tsplot([y, Cctc(:, 1)], datei, names);
                    pause
                    close all
                    if lam == 0
                        names = char('original series in levels ', [compn{i}, ...
                            ' in the original scale']);
                        tsplot([yor, oCctc(:, 1)], datei, names);
                        pause
                        close all
                    end
                elseif ixc > 0 || (i == 5 && irreg == 1)
                    names = compn{i};
                    tsplot(Cctc(:, cont), datei, names);
                    pause
                    close all
                    if lam == 0
                        names = [compn{i}, ' in the original scale'];
                        tsplot(oCctc(:, cont), datei, names);
                        pause
                        close all
                    end
                end
            end
        end
    end
end
