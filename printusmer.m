function printusmer(fid, datei, tname, yor, y, ny, lam, modescr, result, nreg, nbeta)
%*************************************************************************
% This function prints the estimation results of a univariate structural
% model
%
%    INPUTS:
%      fid : file identifier, needed for writing the output into text file
%    datei : calendar structure
%    tname : name of the series (string variable)
%      yor : original time series
%        y : time series used in the estimation etc.
%       ny : length of y
%      lam = 0 : compute logs of y
%          = 1 : do not compute logs
%  modescr : structure containing model information (output of suusm)
%   result : structure with the estimation results (output of usmestim)
%     nreg : number of regression variables in the observation equation;
%    nbeta : number of regression coefficients in the state space model;
%            number of columns of the matrices X and W
%
% Note: nreg corresponds to the number of nonzero columns of X and nbeta to
%       the number of the columns of X and W, where
%       X is an (n x nbeta) matrix containing the X_t matrices;
%         an(1 x nbeta) matrix if it is time invariant; it can be []
%       W is an (n*nalpha x nbeta) matrix containing the W_t matrices;
%         an (nalpha x nbeta) matrix if it is time invariant; it can be []
%       and
%       X_t and W_t are matrices of the state space model:
%
%           y_t = X_t*beta + Z_t*alpha_t + G_t*epsilon_t
%           alpha_{t+1}= W_t*beta + T_t*alpha_t + H_t*epsilon_t
%
%           where epsilon_t is (0,sigma^2I),
%
%           with initial state
%
%           alpha_1= c + W_0*beta + a_1 + A_1*delta
%
%           where c is (0,Omega) and delta is (0,kI) (diffuse)
%
% Copyright (c) January 2010 by Victor Gomez
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
%**************************************************************************

arp = modescr.arp;
cycle = modescr.cycle;
pfix = modescr.pfix;
pvar = modescr.pvar;
if modescr.trend == -1, flevel = 1;
else flevel = 0;
end
if modescr.slope == -1, fslope = 1;
else fslope = 0;
end
if modescr.seas == -1, fseas = 1;
else fseas = 0;
end
xx = modescr.x;
xx(pvar) = result.xvf; %estimated parameters
tt = result.tv;
se = result.se;
hb = result.h;
tb = result.tvr;
seb = result.ser;
sconp = sqrt(result.sigma2c);

%print estimation results
inft = minft(fid, 1, 9, 3, 1);
dateib = datei;
inftb = inft;
fnameb = tname;
prtser(fid, fnameb, yor, y, ny, dateib, inftb, lam);

fprintf(fid, 'Estimation results:\n');
stord = modescr.stord;
conc = modescr.conc;
arpc = arp;
if cycle > 0
    arpc = arp + 2;
end
nst = length(stord) - arpc;
if nst > 0, xx(1:end-arpc) = xx(1:end-arpc) * sconp;
end
sse = zeros(size(xx));
sse(pfix) = NaN;
sse(pvar) = se;
%this line added 12-2011
tte = zeros(size(xx));
tte(pfix) = NaN;
tte(pvar) = tt;
%this line changed 12-2011
xx = [sconp, xx];
sse = [NaN, sse];
tte = [NaN, tte];
z = [xx', sse', tte'];
clear in
in.cnames = char('  Estimate', 'Std. Error', '   T-ratio');
in.rnames = ['Parameter   '];
rnamess = char('Sigma irreg.', 'Sigma level ', 'Sigma slope ', ...
    'Sigma seaso.', 'Sigma autor.', 'Sigma cycle');
in.rnames = [in.rnames; rnamess(conc, :)];
if nst > 0
    for i = 1:nst, in.rnames = [in.rnames; rnamess(stord(i), :)];
    end
end
if cycle > 0
    in.rnames = char(in.rnames, 'Cycle rho ', 'Cycle freq.');
end
if arp > 0
    rnamesa = ['Ar(', num2str(1), ')       '];
    for i = 2:arp, rnamesa = char(rnamesa, ['Ar(', num2str(i), ')       ']);
    end
    for i = 1:arp, in.rnames = [in.rnames; rnamesa(i, :)];
    end
end
in.fmt = char('%12.4f', '%12.4f', '%12.4f');
in.fid = fid;
mprint(z, in);
fprintf(fid, ['Parameter ', rnamess(conc, :), ' is concentrated out of the likelihood\n']);


%print correlations of the estimates
fprintf(fid, '\nCorrelations of the estimates:\n');
idx = ~isnan([NaN, sse]);
rnames = [in.rnames(1, :); in.rnames(idx, :)];
clear in
SS = result.Cv;
DD = diag(sqrt(abs(diag(SS))));
CC = (DD \ SS) / DD;
in.cnames = rnames(2:end, :);
in.rnames = rnames;
in.fmt = '%12.4f';
in.fid = fid;
mprint(CC, in);
fprintf(fid, '\n');

% Coefficient of determination

seas = modescr.seas;
freq = modescr.freq;
% Create vector with all parameters after estimation
xvf = result.xvf;
xf = result.xf;
lx = length(xvf) + length(xf);
x = zeros(lx, 1);
x(pfix) = xf;
x(pvar) = xvf;

sd = zeros(1, 2);
ne = length(result.e); % lenght of the residual vector
Pevf = result.Pevf; % prediction error variance

for i = 2:3
    if any(stord == i)
        j = (stord == i);
        if x(j) == 0 % the estimated or fixed variance is equal to zero
            sd(i-1) = 0;
        else sd(i-1) = 1;
        end
    elseif (conc ~= i) && ~any(stord == i) % the corresponding component is constant
        sd(i-1) = 0;
    else sd(i-1) = 1;
    end
end


if (sd(1) == 0 && sd(2) == 0 && seas == 0)
    rname = 'R^2'; % Data without trend movements
    igood = ~isnan(y);
    imiss = isnan(y);
    gy = y(igood);
    my = mean(gy);
    y(imiss) = my; % Replace missing values of y with the mean of y
    sdy = sum((y - my).^2);
    fact = 1;
elseif ((sd(1) ~= 0 || sd(2) ~= 0) && seas == 0)
    rname = 'R_D^2'; % Data with trend movements
    dy = diff(y);
    igood = ~isnan(dy);
    imiss = isnan(dy);
    gdy = dy(igood);
    mdy = mean(gdy);
    dy(imiss) = mdy; % Replace missing values of dy with the mean of dy
    sdy = sum((dy - mdy).^2);
    % Distinction between the cases without and with explanatory variables
    if nbeta == 0
        fact = 1;
    elseif nbeta > 0
        fact = (ny - 2) / (ne - nreg);
    end
elseif seas ~= 0
    if sd(1) == 0 && sd(2) == 0 % Seasonal data without trend or with constant level
        rname = 'R_S^2';
        % Seasonal mean of y
        smeanf = zeros(freq, 1);
        ly = length(y);
        smean = zeros(ly, 1);
        for i = 1:freq
            sind = i:freq:ly;
            yseas = y(sind); % Data for season i
            igood = ~isnan(yseas);
            imiss = isnan(yseas);
            yseas = yseas(igood);
            smeanf(i) = mean(yseas);
            y(sind(imiss)) = smeanf(i); % Replace missing values of dy in season i with the mean of season i
            smean(sind) = smeanf(i);
        end
        sdy = sum((y - smean).^2);
        % Distinction between the cases without and with explanatory variables
        if nbeta == 0
            fact = 1;
        elseif nbeta > 0
            fact = (ny - freq - 1) / (ne - nreg);
        end
    elseif sd(1) ~= 0 || sd(2) ~= 0
        rname = 'R_S^2'; % Seasonal data with trend
        dy = diff(y);
        % Seasonal mean of dy
        smeanf = zeros(freq, 1);
        ldy = length(dy);
        smean = zeros(ldy, 1);
        for i = 1:freq
            sind = i:freq:ldy;
            dyseas = dy(sind); % Data for season i
            igood = ~isnan(dyseas);
            imiss = isnan(dyseas);
            dyseas = dyseas(igood);
            smeanf(i) = mean(dyseas);
            dy(sind(imiss)) = smeanf(i); % Replace missing values of dy in season i with the mean of season i
            smean(sind) = smeanf(i);
        end
        sdy = sum((dy - smean).^2);
        % Distinction between the cases without and with explanatory variables
        if nbeta == 0
            fact = 1;
        elseif nbeta > 0
            fact = (ny - freq - 1) / (ne - nreg);
        end
    end
end

% Formulas for R^2, Harvey(1989) "Forecasting, structural time series
% models and the Kalman filter":
% eq.(5.5.14), eq.(5.5.17) for models without explanatory variables
% eq.(7.4.11), eq.(7.4.12) for models with explanatory variables,

Rsq = 1 - (ne * Pevf * fact) / sdy;
% *************************************************************************
fprintf(fid, ['\nCoefficient of determination ', rname, ' %11.4f'], Rsq);
fprintf(fid, '\nPrediction error variance          %11.4f', Pevf);

SPevf = result.SPevf;
fprintf(fid, '\nResidual standard error            %11.4f\n', SPevf);

%print regression variables
if nbeta > 0
    fprintf(fid, '\n');
    fprintf(fid, 'Regression parameters:\n');
    %  seb=sqrt(diag(Mb*conp)); tb=hb./seb;             %standard errors and t-values
    Mbeta = [hb, seb, tb];
    clear in
    in.cnames = char('  Estimate', 'Std. Error', '   T-ratio');
    rnames = char('Parameter');
    cont = 0;
    if flevel == 1
        rnames = char(rnames, ['level']);
        cont = cont + 1;
    end
    if fslope == 1
        rnames = char(rnames, ['slope']);
        cont = cont + 1;
    end
    if fseas == 1
        for i = 1:freq - 1
            rnames = char(rnames, ['seas', num2str(i)]);
            cont = cont + 1;
        end
    end;
    mreg = nbeta - cont;
    if mreg > 0
        for i = 1:mreg
            rnames = char(rnames, ['reg', num2str(i)]);
        end
    end
    in.rnames = rnames;
    in.fmt = char('%12.5f', '%12.5f', '%8.2f');
    in.fid = fid;
    mprint(Mbeta, in);
    fprintf(fid, '\n');
    
    %print correlations of the regression estimates
    fprintf(fid, '\nCorrelations of the regression estimates:\n');
    clear in
    M = result.M;
    DDr = diag(sqrt(abs(diag(M))));
    CCr = (DDr \ M) / DDr;
    in.cnames = rnames(2:end, :);
    in.rnames = rnames;
    in.fmt = '%12.4f';
    in.fid = fid;
    mprint(CCr, in);
    fprintf(fid, '\n');
    
end

if isfield(modescr, 'nmiss')
    nmiss = modescr.nmiss;
    if nmiss > 0
        Interp = modescr.Interp;
        sInterp = modescr.sInterp;
        idxn = modescr.idxn;
        clear in
        in.cnames = char('  Estimate', 'Std. Error');
        rnamesrgi = ['interp. ', num2str(idxn(1))];
        for i = 2:nmiss
            rnamesrgi = char(rnamesrgi, ['interp. ', num2str(idxn(i))]);
        end
        rnames = char('Interpolated value', rnamesrgi);
        in.rnames = rnames;
        in.fmt = char('%12.5f');
        in.fid = fid;
        mprint([Interp, sInterp], in);
        fprintf(fid, '\n');
        if (lam == 0)
            oInterp = modescr.oInterp;
            osInterp = modescr.osInterp;
            rnames = char('Interpolated value (levels)', rnamesrgi);
            in.rnames = rnames;
            mprint([oInterp, osInterp], in);
            fprintf(fid, '\n');
        end
    end
end
