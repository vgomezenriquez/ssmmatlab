%script file containing ARIMA default values

if isfield(ser, 'Y') %Y contains the regression variables
    Y = ser.Y;
    if isfield(ser, 'rnamesrg') %names for regression variables
        rnamesrg = ser.rnamesrg;
    else
        rnamesrg = [];
    end
    if ~isempty(rnamesrg)
        rnamesr = 1;
    else
        rnamesr = 0;
    end
else
    Y = [];
    rnamesr = 0;
    rnamesrg = [];
end

if isfield(ser, 'npr')
    npr = ser.npr; %number of forecasts
else
    npr = 0;
end

if isfield(ser, 'pr') %printing results in external file
    pr = ser.pr; % pr = 1 print in external file
else %    = 0 do not print in external file
    pr = 1;
end


%transfer function
if isfield(ser, 'ninput')
    ninput = ser.ninput; %ninput is the number of input variables
    if isfield(ser, 'Yin') %Yin contains the input variables
        Yin = ser.Yin;
    else
        error('there must be a field Yin for the input variables')
    end
    if isfield(ser, 'prelivar')
        prelivar = ser.prelivar;
    else
        prelivar = 0;
    end
    if isfield(ser, 'tfident')
        tfident = ser.tfident; %tfident is the parameter for automatic identification
    else
        tfident = 0; %default value
    end
    if (tfident == 0 && prelivar == 0)
        if ~isfield(ser, 'delay') || ~isfield(ser, 'ma') || ~isfield(ser, 'ar')
            error('fields delay, ma and ar should be input if tfident = 0')
        end
    end
    %a negative nlagtf means that the program will compute the
    %number of lags. A certain number of lags in the inputs is needed to
    %perform automatic model identification, outlier detection, etc. in the
    %first run of the program.
    if isfield(ser, 'nlagtf')
        nlagtf = ser.nlagtf;
    else
        nlagtf = -1;
    end
    if isfield(ser, 'inc')
        inc = ser.inc;
    else
        inc = 1;
    end
    if isfield(ser, 'rnamesi') %rnamesi=1, names for the input variables
        rnamesi = ser.rnamesi;
    else
        rnamesi = 0;
    end
    if isfield(ser, 'rnamesiv') %names for the input variables
        rnamesiv = ser.rnamesiv;
    else
        rnamesiv = [];
    end
    if (npr > 0)
        if ~isfield(ser, 'modinput') || ~isfield(ser, 'modpred')
            error('If npr > 0, structure fields modinput and modpred should be input')
        end
        mYin = size(Yin, 2);
        for i = 1:mYin
            if ~isfield(ser.modinput(i), 'mod')
                error('If npr > 0, field modinput(i).mod should be input')
            elseif ser.modinput(i).mod == 1 && (~isfield(ser.modinput(i), 'phi') ...
                    || ~isfield(ser.modinput(i), 'theta') || ~isfield(ser.modinput(i), 'sigma2'))
                error('If modinput(i).mod is one, fields phi, theta and sigma2 should be input')
            end
            if ~isfield(ser.modpred(i), 'pred')
                error('If npr > 0, field modpred(i).pred should be input')
            elseif size(ser.modpred(i).pred, 2) ~= 1
                error('If npr > 0, field modpred(i).pred should be an (npr x 1) array')
            end
        end
    end
else
    ninput = 0;
    nlagtf = 0;
    tfident = 0;
    Yin = [];
    inc = 1;
    rnamesi = 0;
    rnamesiv = [];
end
[nYin, mYin] = size(Yin);
ycii = [];
%transfer function

ny = length(yor); %series length
[nY, mY] = size(Y);
nreg = mY; %nreg = number of initial regression parameters.
%The mean is not included

%generation of names for regression variables
if (rnamesr == 0) && (nreg > 0)
    rnamesrg = ['reg', num2str(1)];
    for i = 2:nreg
        rnamesrg = char(rnamesrg, ['reg', num2str(i)]);
    end
    rnamesr = 1;
end

%transfer function
%checking for inconsistencies
if (ninput ~= mYin)
    error('ninput should be equal to the number of columns in Yin')
end
%end of transfer function

% yin=[]; % we will store here the original series when there are lagged variables
% nyin=0;
Yreg = Y; % and here the original regression variables
initreg = nreg;
chb = 0; %chb=1, in fstlkhev OLS estimator and MSE are computed

%transfer function
myl = 0;
if ninput > 0
    if prelivar == 1
        %
        %preliminary VAR analysis
        %
        prtx = 1;
        minlags = 0;
        a = 1.6;
        pt = ser.freq + 5;
        maxlags = ceil(max(log(ny)^a, pt));
        lagsopt = lratiocr([yor, Yin], maxlags, minlags, prtx);
        flagsopt = 0;
        if lagsopt == 0
            lagsopt = 1;
            flagsopt = 1;
        end
        %estimate var
        resv = var_est([yor, Yin], lagsopt, 1);
        outa.varlag = lagsopt;
        outa.varstr = resv;
        if (pr >= 1)
            disp(' ')
            if flagsopt == 0
                disp('Estimated order in VAR:  ')
                disp(lagsopt)
            else
                disp('Estimated order in VAR is zero')
                disp('VAR order changed to one')
            end
            disp(' ');
            disp('***** Estimated VAR Model  *****');
            disp(' ');
            clear in
            in.fid = 1;
            in.fmt = char('%12.4f');
            tit = 'AR';
            strt = 1;
            mprintar(resv.phi(:, :, 2:lagsopt+1), in, tit, strt);
            disp(' ')
            tit = 'Constant';
            mprintar(resv.const', in, tit);
            
            disp(' ');
            disp('***** Estimated t-values  *****');
            disp(' ');
            clear in
            in.fid = 1;
            in.fmt = char('%12.4f');
            tit = 'tv-AR';
            strt = 1;
            mprintar(resv.phitv(:, :, 2:lagsopt+1), in, tit, strt);
            disp(' ')
            tit = 'tv-Constant';
            mprintar(resv.consttv', in, tit);
            SS = resv.sigmar;
            disp(' ');
            tit = 'Sigma:';
            mprintar(SS, in, tit);
            DD = diag(sqrt(abs(diag(SS))));
            CC = (DD \ SS) / DD;
            disp(' ');
            tit = 'Correlations:';
            mprintar(CC, in, tit);
        end
        %
        %end of preliminary VAR analysis
        %
        return
    end
    yin = yor;
    chb = 1;
    nyin = length(yin);
    [nYin, mYin] = size(Yin);
    %  nreg=nreg+ninput;
    if (nlagtf < 0) %we use a formula for the number of lags in the inputs
        nlagtf = min(13, ceil(log(nyin)^1.5));
    end
    rnamesr = 1; %flag for regression variables names
    if rnamesi == 1 %names for input variables given
        for j = 1:ninput
            for i = 0:nlagtf
                rnamesrg = strvcat(rnamesrg, [rnamesiv(j, :), '_', num2str(i)]);
            end
        end
    else
        for j = 1:ninput %no names given for input variables
            for i = 0:nlagtf
                rnamesrg = strvcat(rnamesrg, ['inp', num2str(j), '_', num2str(i)]);
            end
        end
    end
    ylags = []; %set up input lags
    for i = 1:ninput
        xlags = ltflag(Yin(:, i), nlagtf);
        ylags = [ylags, xlags];
    end
    mninput = nlagtf + 1; %initial observations lost due to lags
    [nyl, myl] = size(ylags);
    yor = yor(mninput:end);
    %eliminate dummy variables corresponding to missing observations in the
    %first observations (before nlagtf+1)
    npmiss = 0;
    rnamesrgm0 = rnamesrgm;
    if nmiss > 0
        present = find(idxn < mninput);
        npmiss = length(present);
        if (npmiss > 0)
            npmiss1 = npmiss + 1;
            if npmiss1 > mY
                Y = [];
                rnamesrgm = [];
            else
                Y = Y(:, npmiss1:mY);
                rnamesrgm = rnamesrgm(npmiss1:end, :);
            end
            rnamesrg = rnamesrg(npmiss1:end, :);
            mY = mY - npmiss;
            initreg = initreg - npmiss;
        end
    end
    if (mY > 0)
        Y = [Y(mninput:end-npr, 1:mY), ylags]; % put the lags after the regression variables
    else
        Y = ylags;
    end
    datei = cal(bg_year, bg_per, freq);
    res = cal(datei, mninput); % compute the new initial date due to the input lags
    bg_year = res.year;
    bg_per = res.period;
    ny = length(yor); %series length
    [nY, mY] = size(Y);
    nreg = mY;
    %nreg = number of initial regression parameters, included the lags
    nelinput = ones(ninput, 1) * (nlagtf + 1); %number of weights for each input variable
    cnelinput = nelinput;
    ninput = 0;
end
%end of transfer function

%series default parameters
s = freq;
%Arima orders
if isfield(ser, 'dr')
    dr = ser.dr;
else
    dr = 1;
end
if isfield(ser, 'ds')
    ds = ser.ds;
else
    ds = 1;
end
if isfield(ser, 'p')
    p = ser.p;
else
    p = 0;
end
if isfield(ser, 'ps')
    ps = ser.ps;
else
    ps = 0;
end
if isfield(ser, 'q')
    q = ser.q;
else
    q = 1;
end
if isfield(ser, 'qs')
    qs = ser.qs;
else
    qs = 1;
end
if isfield(ser, 'S')
    S = ser.S;
    if isfield(ser, 'dS')
        dS = ser.dS;
    else
        dS = 1;
    end
    if isfield(ser, 'qS')
        qS = ser.qS;
    else
        qS = 1;
    end
else
    S = 0;
    dS = 0;
    qS = 0;
end
if s <= 1, ds = 0;
    ps = 0;
    qs = 0;
    S = 0;
    dS = 0;
    qS = 0;
end


%  flagm=0;                   % mean for the series (=1, put mean, =0 no mean)
if isfield(ser, 'flagm')
    flagm = ser.flagm;
else
    flagm = 0;
end

nr = p + ps + q + qs + qS; %number of arma parameters
fvdif = 0;
c1 = 0;
maxr = 2;
maxpq = 3; %max regular diff. max regular p and q.
%  autmid=1;                  %automatic model identification of ARIMA model
if isfield(ser, 'autmid')
    autmid = ser.autmid;
else
    autmid = 1;
end
%  fixdif=0;                  %differencing orders fixed (=1)
if isfield(ser, 'fixdif')
    fixdif = ser.fixdif;
else
    fixdif = 0;
end
%  pfix=[]; pvar=1:nr;        %fixed and free parameters
%  xi=zeros(1,nr);
%
%Fixing of parameters: There should be two fields, pfix and vfix, in
%the ser structure. The first one for the parameter indices and the second
%one for the parameter values to be fixed.
%The order of the parameters should be: 1,..p,p+1,..p+ps,p+ps+1,..p+ps+q,
%p+ps+q+1,..,p+ps+q+qs,p+ps+q+qs+1,p+ps+q+qs+qS.
%The order in pfix should coincide with that in vfix
%
if isfield(ser, 'pfix')
    pfix = ser.pfix;
    if isfield(ser, 'vfix')
        %vfix should be a row array containing the fixed parameter values
        vfix = ser.vfix;
    else
        error('if parameters are fixed, vfix should be a field')
    end
    tidx = 1:nr;
    tidx(pfix) = [];
    pvar = tidx;
    xi = zeros(1:nr);
    xi(pfix) = vfix;
    if autmid == 1
        disp('if parameters are fixed, no automatic model identification can be performed')
        error('set ser.autmid to zero when fixing parameters')
    end
else
    pfix = [];
    pvar = 1:nr;
    xi = zeros(1, nr);
end
%  lam=1;                    %logs of data (=0 logs are taken, <0 automatically identified)
if isfield(ser, 'lam')
    lam = ser.lam;
else
    lam = -1;
end

fmarqdt = 1; %parameter for estimation method (1=marqrt)
fjac = 0; %parameter for method to obtain standar error (1=jacobian)
cw = 1.96; %coefficient for the confidence bands
ols = 0;
a = 2.0; %1.5;         %parameters for the Hannan-Rissanen method
% ols = 0 use the Levinson-Durbin algorithm
%     = 1 use OLS
% TRAMO uses ols = 0
% a = the exponent in log(n)^a for the length of the long AR
% TRAMO uses a=2.0

% delta = value for TC outliers
delta = .7;
%  out=1;                   % out = 1 perform outlier detection
%     = 0 do not perform outlier detection
if isfield(ser, 'out')
    out = ser.out;
else
    out = 0;
end
%  omet=1;                    % omet = 1 use exact ML for model estimation
%      = 0 use Hannan-Rissanen
if isfield(ser, 'omet')
    omet = ser.omet;
else
    omet = 0;
end
% C    = critical value for outlier detection
%      if negative, it is computed
%      depending on the sample size
if isfield(ser, 'C')
    C = ser.C;
else
    C = -1;
end
%default values of C0. Also the values for automatic model identification
%if ny <= 50, C0=2.8; elseif (ny> 50 & ny <= 75), C0=3.2; elseif (ny> 75 & ny <= 100),...
%C0=3.5; elseif (ny > 100 & ny <= 200), C0=3.8; elseif (ny >200 & ny<=300), C0=3.9;...
%elseif (ny >300 & ny<=400), C0=4.; elseif (ny >400 & ny<=500), C0=4.1; else, C0=4.2; end
%C0 is used in the log test and automatic model identification
if isfield(ser, 'C0')
    C0 = ser.C0;
else
    C0 = 2.6 + log(log(ny));
end
% schr = 0 outliers of type AO and TC are considered (default)
%      = 1 outliers of type AO, TC and LS are considered
if isfield(ser, 'schr')
    schr = ser.schr;
else
    schr = 1;
end
%  (sp1,sp2)                    %span for outlier detection
if isfield(ser, 'sp1')
    sp1 = ser.sp1;
else
    sp1 = 1;
end
if isfield(ser, 'sp2')
    sp2 = ser.sp2;
else
    sp2 = ny;
end


%  trad=0;                   %trading day (TD) effect (trad=-1, test)
%  tradval=[1 6];            %tradval=possible number of TD variables  (0 is also a value)
%  leapy=0;                  %leap year (leapy=-1, test)
%  easte=0;                  %easter effect (easte=-1, test)
%  durval=[4 6];             %durval=possible days previous to Easter (0 is also a value)
if isfield(ser, 'trad')
    trad = ser.trad;
else
    trad = 0;
end
if isfield(ser, 'tradval')
    tradval = ser.tradval;
else
    tradval = [1, 6];
end
if isfield(ser, 'leapy')
    leapy = ser.leapy;
else
    leapy = 0;
end
if isfield(ser, 'easte')
    easte = ser.easte;
else
    easte = 0;
end
if isfield(ser, 'durval')
    durval = ser.durval;
else
    durval = [4, 6];
end


lag = min(36, max([floor(.15*ny), ...
    3 * s, 10]));%number of lags for autocorrelations

if isfield(ser, 'nlestim')
    nlestim = ser.nlestim; %nonlinear estimation
else
    nlestim = 1;
end

if isfield(ser, 'mvx') %mvx :  =1 exact maximum likelihood
    mvx = ser.mvx; %       =0 unconditional least squares
else
    mvx = 1;
end

if isfield(ser, 'olsres') %olsres :  =1 use OLS residuals
    olsres = ser.olsres; %          =0 do not use OLS residuals
else
    olsres = 0;
end


f = 'fasttf';
tr = 0; %parameters for Levenberg-Marquardt method:
tolf = 1e-4; %f  :   a function to evaluate the vector ff of individual
maxit = 100;
nu0 = .0; %       functions such that ff'*ff is minimized
jac = 0;
prt = 0; %tr :   >0 x is passed from marqdt to f but not passed from f
%       to marqdt
%       =0 x is passed from marqdt to f and passed from f to
%       marqdt
%tolf:  a parameter used for stopping
%       TRAMO uses tolf=1e-4
%jac:   =1 evaluation of jacobian at gradient at the solution
%       is performed
%       =0 no evaluation of jacobian at gradient at the
%       solution is performed
%maxit: maximum number of iterations
%nu0:   initial value of the nu parameter
%prt:   =1 printing of estimation details
%       =2 printing of more estimation details
%       =0 no printing of estimation details


fh = 1;
wd = 9;
nd = 2;
scale = 1; %table parameters:
%fh: flag for header and years
%wd: format width
%nd: number of decimal points
%scale: =1 scale data if necessary
%       =0 do not scale data

%  gft=1;                      %flag for plots:
% gft = 1 plot series at the end
%     = 0 do not plot series at the end
if isfield(ser, 'gft')
    gft = ser.gft;
else
    gft = 0;
end

gf11 = 1; %flag for one by one plots:
% gf11 = 1 pause after each plot
%      = 0 do not pause after each plot

%parameters for transfer function
if isfield(ser, 'maxndtf') %maximum degree for numerator in transfer function indentification
    maxndtf = ser.maxndtf;
else
    maxndtf = 2;
end
if isfield(ser, 'maxddtf') %maximum degree for denominator in transfer function identification
    maxddtf = ser.maxddtf;
else
    maxddtf = 2;
end
if isfield(ser, 'backwd') %backward elimination for transfer function identification
    backwd = ser.backwd;
else
    backwd = 1;
end
if isfield(ser, 'Cb') %critical value for backward elimination in transfer function identification
    Cb = ser.Cb;
else
    Cb = 2.;
end

%after generating the program figures, note that:
%1) a call to close all deletes all graphic figures
%2) delete(1:5) deletes figures 1 to 5
%3) figure(3) makes visible figure 3
