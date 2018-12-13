function [out, ser] = tfeasy(y, x, freq, varargin)
%**********************************************************************
%                       EASY TRANSFER FUNCTION MODELING
%
%                            USAGE :
% out = tfeasy(y,x,freq,'option1',optionvalue1,'option2',optionvalue2,...)
%
%       INPUTS :
%------------------
%      REQUIRED
%            y : (ly x 1) array containing the series;
%            x : (ly x ni) array containing the inputs
%         freq : data frequency (number of observations per year)
%------------------
%       OPTIONS
%  '[bg_year bg_per]': (1 x 2) array containing the initial year and the
%                       initial period. Default [2000 1]
%               'lam': data transformation (logs), = 0 logs, =1 no logs,
%                      default -1 (test for logs)
%          '[p dr q]': (1 x 3) array containing the regular orders
%                      default: [0 1 1]
%        '[ps ds qs]': (1 x 3) array containing the first seasonal orders
%                      default: [0 1 1]
%                 'S': second seasonality. Default 0
%           '[dS qS]': (1 x 2) array containing the second seasonal orders
%                      default: [1 1]
%             'flagm': flag for mean, =1 mean, =0, no mean, default 0
%                      It has not effect with automatic model
%                      identification
%              'pfix': index array for fixed parameters
%              'vfix': array for fixed parameter values
%            'fixdif': flag for fixing the differencing degrees, =1
%                      degrees are fixed, = 0 not fixed, default 0
%            'autmid': flag for automatic model identification, = 1,
%                      perform automatic model identification, = 0, no
%                      automatic model identification, default 1
%                 'Y': array for regression variables, default []
%          'rnamesrg': string matrix for names of regression variables,
%                      default []
%           'nlestim': flag for nonlinear estimation, = 1, nl estimation,
%                      = 0, no nl estimation, default 1
%               'mvx': flag for nl method, = 1, exact maximum likelihood,
%                      = 0, unconditional least squares, default 1
%               'npr': number of forecasts, default 0
%            'olsres': flag for OLS residuals, = 1, OLS residuals are used,
%                      = 0, uncorrelated residuals (transformation of OLS
%                      residuals) are used, default 0
%                'pr': flag for printing in an external file, = 1, printing
%                       = 0, no printing, default 1
%               'gft': flag for graphics, = 1, plot series, = 0, no plots
%                       = 2, plots are saved but not displayed, = 3, plots
%                       are both saved and displayed, default 0
%               'out': out = 1 perform outlier detection
%                      = 0 do not perform outlier de
%              'omet': omet = 1 use exact ML for model estimation
%                      = 0 use Hannan-Rissanen
%                 'C': critical value for outlier detection
%                      if negative, it is computed depending on the
%                      sample size
%                'C0': critical value for outlier detection used in the log
%                      test and automatic model identification, default
%                      C0=2.6 + log(log(ny)) (ny = series length)
%              'schr': = 0 outliers of type AO and TC are considered, =1
%                      outliers of type AO, TC and LS are considered,
%                      default 1
%               'sp1': (sp1,sp2) span for outlier detection, default sp1 =1
%                      default sp2=ny, where ny = series length
%               'sp2':
%              'trad': = 0 no trading day effect, = 1 TD effect, = -1, test
%                      for TD effect, default 0
%           'tradval': possible number of TD variables (0 is also a value),
%                      default [1 6]
%             'leapy': = 0, no leap year effect, = 1 LP effect, = -1, test
%                      for LP effect, default 0
%             'easte': = 0 no Easter effect, = 1 Easter effect, = -1, test
%                      for Easter effect, default 0
%            'durval': possible days previous to Easter (0 is also a value)
%                      default [4 6]
%             'sname': character array containing the series name
%                      default series1
%           'rnamesi': flag for names of the input variables, = 1, names are
%                      given by the user, =0, names are given by the program,
%                      default 0
%          'rnamesiv': character array containing the names for the input
%                      variables, default []
%          'prelivar': flag for preliminary VAR analysis, = 1, perform VAR
%                      analysis, = 0, no VAR analysis, default 0
%             'delay': array containing the filter delays, if tfident=0 and
%                      prelivar=0
%                'ma': array containing the filter ma degrees, if tfident=0
%                      and prelivar=0
%                'ar': array containing the filter ar degrees, if tfident=0
%                      and prelivar=0
%               'inc': = 0, the initial states in the filter equations to obtain
%                      the filtered variables are equal to zero (not estimated)
%                      = 1, the initial states in the filter equations are
%                      estimated
%          'modinput': structure containing the input models in subfields
%                      phi, theta and sigma2 if subfield mod = 1; default
%                      mod 0. The input model is used to compute the mse,
%                      not the input forecasts. It should contain the
%                      nonstationary part.
%           'modpred': structure containing the input forecasts in subfield
%                      pred. If npr > 0, the user should provide the input
%                      forecasts as an (npr x 1) array for each input,
%                      whether there is a model for the input or not
%           'tfident': flag for automatic TF identification, default 0
%            'backwd': flag for backward elimination in transfer function
%                      identification, default 0
%                'Cb': critical value for backwar elimination in transfer
%                      function identification, default 2.
%            'nlagtf': number of lags for automatic model identification. If
%                      negative, the program will compute the number of lags.
%                      default, -1
%           'maxndtf': maximum degree for numerator in transfer function
%                      identification
%           'maxddtf': maximum degree for denominator in transfer function
%                      identification%
%
%       OUTPUT : a structure, the output of function arimaestni
%------------------
%
%     Examples:
%
%   [out,ser]=tfeasy(y,x,freq,'tfident',1,'out',1)
%   out=arimaeasy(y,x,freq,'[p dr q]',[0 1 1],'leapy',-1,'tfident',1)
%
%
% Copyright (c) February 2018 by Victor Gomez
% Ministerio de Hacienda y A.P., Direccion Gral. de Presupuestos,
% Subdireccion Gral. de Analisis y P.E.,
% Alberto Alcocer 2, 1-P, D-34, 28046, Madrid, SPAIN.
% Phone : +34-915835439
% E-mail: vgomez@sepg.minhap.es
% The author assumes no responsibility for errors or damage resulting from the use
% of this code. Usage of this code in applications and/or alterations of it should
% be referenced. This code may be redistributed if nothing has been added or
% removed and no money is charged. Positive or negative feedback would be appreciated.
%
%*************************************************************************

% Check the number and validity of arguments and set defaults

if nargin < 3
    error('There must be at least three inputs to tfeasy');
end
[ny, my] = size(y);
if my > 1
    error('Array y must be (n x 1) in arimaeasy');
end
Yin = x;
[nYin, mYin] = size(Yin);
if nYin ~= ny
    error('Input and output must have the same numer of observations in tfeasy')
end
if (~isscalar(freq) || (freq < 0))
    error('Value of freq must be a positive integer in arimaeasy');
end


ser.yor = y;
ser.Yin = Yin;
ser.ninput = mYin;
ser.freq = freq;

% set default values for '[bg_year bg_per]'
def = cell(1, 47);
ldef = length(def);
for i = 1:ldef
    def{i} = [];
end
def{1} = [2000, 1];
ser.bg_year = 2000;
ser.bg_per = 1;
% set default values
def{32} = 'mtfseries';
sname = 'mtfseries';

if nargin == 3
    val = def;
end

lvar = length(varargin);

if lvar > 0
    if mod(lvar, 2) == 1 % number of elements is odd
        error('To specify an option, one of the admissible types and an admissible value must be chosen');
    else
        opt = {'[bg_year bg_per]', ... %1
            'lam', ... %2
            '[p dr q]', ... .          %3
            '[ps ds qs]', ... %4
            'S', ... %5
            '[dS qS]', ... %6
            'flagm', ... %7
            'pfix', ... %8
            'vfix', ... %9
            'fixdif', ... %10
            'autmid', ... %11
            'Y', ... %12
            'rnamesrg', ... %13
            'nlestim', ... %14
            'mvx', ... %15
            'npr', ... %16
            'olsres', ... %17
            'pr', ... %18
            'gft', ... %19
            'out', ... %20
            'omet', ... %21
            'C', ... %22
            'C0', ... %23
            'schr', ... %24
            'sp1', ... %25
            'sp2', ... %26
            'trad', ... %27
            'tradval', ... %28
            'leapy', ... %29
            'easte', ... %30
            'durval', ... %31
            'sname', ... %32
            'rnamesi', ... %33
            'rnamesiv', ... %34
            'prelivar', ... %35
            'delay', ... %36
            'ma', ... %37
            'ar', ... %38
            'inc', ... %39
            'modinput', ... %40
            'modpred', ... %41
            'tfident', ... %42
            'backwd', ... %43
            'Cb', ... %44
            'nlagtf', ... %45
            'maxndtf', ... %46
            'maxddtf'};%47
        ldef = length(opt);
        val = def;
        for i = 1:2:(lvar - 1)
            if any(strcmpi(varargin{i}, opt))
                for k = 1:ldef
                    if strcmpi(varargin{i}, opt{k})
                        val{k} = varargin{i+1};
                    end
                end
            else
                error('Option must be of an admissible type');
            end
        end
    end
end
% Check the validity of option values if they are not set to defaults

for i = 1:ldef
    if ~isequal(val{i}, def{i})
        if strcmpi('[bg_year bg_per]', opt{i})
            if (~isnumeric(val{i}) || (size(val{i}, 1) ~= 1 || size(val{i}, 2) ~= 2))
                error('Value of option [bg_year bg_per] must be a numeric (1 x 2) array');
            else
                ser.bg_year = val{i}(1);
                ser.bg_per = val{i}(2);
            end
        elseif strcmpi('lam', opt{i})
            if (~isscalar(val{i}) || (val{i} ~= 0 && val{i} ~= 1 && val{i} ~= -1))
                error('Value of option lam must be zero, one or minus one');
            else
                ser.lam = val{i};
            end
        elseif strcmpi('[p dr q]', opt{i})
            if (~isnumeric(val{i}) || (size(val{i}, 1) ~= 1 || size(val{i}, 2) ~= 3))
                error('Value of option [p dr q] must be a numeric (1 x 3) array');
            else
                ser.p = val{i}(1);
                ser.dr = val{i}(2);
                ser.q = val{i}(3);
            end
        elseif strcmpi('[ps ds qs]', opt{i})
            if (~isnumeric(val{i}) || (size(val{i}, 1) ~= 1 || size(val{i}, 2) ~= 3))
                error('Value of option [ps ds qs] must be a numeric (1 x 3) array');
            else
                ser.ps = val{i}(1);
                ser.ds = val{i}(2);
                ser.qs = val{i}(3);
            end
        elseif strcmpi('S', opt{i})
            if (~isscalar(val{i}) || (val{i} <= 1 || val{i} == freq))
                error('Value of option S must be greater than one and different from freq');
            else
                ser.S = val{i};
            end
        elseif strcmpi('[dS qS]', opt{i})
            if (~isnumeric(val{i}) || (size(val{i}, 1) ~= 1 || size(val{i}, 2) ~= 2 ...
                    || val{i}(1) > 1 || val{i}(2) > 1))
                error('Value of option [dS qS] must be a numeric (1 x 2) array with elements not greater than one');
            else
                ser.dS = val{i}(1);
                ser.qS = val{i}(2);
            end
        elseif strcmpi('flagm', opt{i})
            if (~isscalar(val{i}) || (val{i} ~= 0 && val{i} ~= 1))
                error('Value of option flagm must be zero or one');
            else
                ser.flagm = val{i};
            end
        elseif strcmpi('pfix', opt{i})
            if ~isnumeric(val{i})
                error('Value of option pfix must be an index array');
            else
                ser.pfix = val{i};
            end
        elseif strcmpi('vfix', opt{i})
            if ~isnumeric(val{i})
                error('Value of option vfix must be numeric array');
            else
                ser.vfix = val{i};
            end
        elseif strcmpi('fixdif', opt{i})
            if (~isscalar(val{i}) || (val{i} ~= 0 && val{i} ~= 1))
                error('Value of option fixdif must be zero or one');
            else
                ser.fixdif = val{i};
            end
        elseif strcmpi('autmid', opt{i})
            if (~isscalar(val{i}) || (val{i} ~= 0 && val{i} ~= 1))
                error('Value of option autmid must be zero or one');
            else
                ser.autmid = val{i};
            end
        elseif strcmpi('Y', opt{i})
            if ~isnumeric(val{i})
                error('Value of option Y must be a numeric array');
            else
                ser.Y = val{i};
            end
        elseif strcmpi('rnamesrg', opt{i})
            if ~ischar(val{i})
                error('Value of option rnamesrg must be a character array');
            else
                if ~isfield(ser, 'Y')
                    error('There must be option Y with option rnamesrg')
                elseif size(val{i}, 1) ~= size(ser.Y, 2)
                    error('Number of rnamesrg must be equal to the number of columns of Y');
                end
                ser.rnamesrg = val{i};
            end
        elseif strcmpi('nlestim', opt{i})
            if (~isscalar(val{i}) || (val{i} ~= 0 && val{i} ~= 1))
                error('Value of option nlestim must be zero or one');
            else
                ser.nlestim = val{i};
            end
        elseif strcmpi('mvx', opt{i})
            if (~isscalar(val{i}) || (val{i} ~= 0 && val{i} ~= 1))
                error('Value of option mvx must be zero or one');
            else
                ser.mvx = val{i};
            end
        elseif strcmpi('npr', opt{i})
            if (~isscalar(val{i}) || (val{i} < 0))
                error('Value of option npr must be zero or positive');
            else
                ser.npr = val{i};
            end
        elseif strcmpi('olsres', opt{i})
            if (~isscalar(val{i}) || (val{i} ~= 0 && val{i} ~= 1))
                error('Value of option olsres must be zero or one');
            else
                ser.olsres = val{i};
            end
        elseif strcmpi('pr', opt{i})
            if (~isscalar(val{i}) || (val{i} ~= 0 && val{i} ~= 1))
                error('Value of option pr must be zero or one');
            else
                ser.pr = val{i};
            end
        elseif strcmpi('gft', opt{i})
            if (~isscalar(val{i}) || (val{i} ~= 0 && val{i} ~= 1 && val{i} ~= 2 && val{i} ~= 3))
                error('Value of option gft must be zero, one, two or three');
            else
                ser.gft = val{i};
            end
        elseif strcmpi('out', opt{i})
            if (~isscalar(val{i}) || (val{i} ~= 0 && val{i} ~= 1))
                error('Value of option out must be zero or one');
            else
                ser.out = val{i};
            end
        elseif strcmpi('omet', opt{i})
            if (~isscalar(val{i}) || (val{i} ~= 0 && val{i} ~= 1))
                error('Value of option omet must be zero or one');
            else
                ser.omet = val{i};
            end
        elseif strcmpi('C', opt{i})
            if ~isscalar(val{i})
                error('Value of option C must be a scalar');
            else
                ser.C = val{i};
            end
        elseif strcmpi('C0', opt{i})
            if ~isscalar(val{i})
                error('Value of option C0 must be a scalar');
            else
                ser.C0 = val{i};
            end
        elseif strcmpi('schr', opt{i})
            if (~isscalar(val{i}) || (val{i} ~= 0 && val{i} ~= 1))
                error('Value of option schr must be zero or one');
            else
                ser.schr = val{i};
            end
        elseif strcmpi('sp1', opt{i})
            if (~isscalar(val{i}) || (val{i} < 0) || (val{i} > ny))
                error('Value of option sp1 must be a positive integer smaller than the series length');
            else
                ser.sp1 = val{i};
            end
        elseif strcmpi('sp2', opt{i})
            if (~isscalar(val{i}) || (val{i} < 0) || (val{i} > ny))
                error('Value of option sp2 must be a positive integer smaller than the series length');
            else
                ser.sp2 = val{i};
            end
        elseif strcmpi('trad', opt{i})
            if (~isscalar(val{i}) || (val{i} ~= 0 && val{i} ~= 1 && val{i} ~= -1))
                error('Value of option trad must be zero, one or minus one');
            else
                ser.trad = val{i};
            end
        elseif strcmpi('tradval', opt{i})
            if (~isnumeric(val{i}) || (size(val{i}, 1) ~= 1))
                error('Value of option tradval must be a numeric array');
            else
                ser.tradval = val{i};
            end
        elseif strcmpi('easte', opt{i})
            if (~isscalar(val{i}) || (val{i} ~= 0 && val{i} ~= 1 && val{i} ~= -1))
                error('Value of option easte must be zero, one or minus one');
            else
                ser.easte = val{i};
            end
        elseif strcmpi('durval', opt{i})
            if (~isnumeric(val{i}) || (size(val{i}, 1) ~= 1))
                error('Value of option durval must be a numeric array');
            else
                ser.durval = val{i};
            end
        elseif strcmpi('leapy', opt{i})
            if (~isscalar(val{i}) || (val{i} ~= 0 && val{i} ~= 1 && val{i} ~= -1))
                error('Value of option leapy must be zero, one or minus one');
            else
                ser.leapy = val{i};
            end
        elseif strcmpi('sname', opt{i})
            if ~ischar(val{i})
                error('Value of option sname must be a character array');
            else
                sname = val{i};
            end
        elseif strcmpi('rnamesi', opt{i})
            if (~isscalar(val{i}) || (val{i} ~= 0 && val{i} ~= 1))
                error('Value of option rnamesi must be zero or one');
            else
                ser.rnamesi = val{i};
            end
        elseif strcmpi('rnamesiv', opt{i})
            if ~ischar(val{i})
                error('Value of option rnamesiv must be a character array');
            else
                if size(val{i}, 1) ~= size(x, 2)
                    error('Number of rnamesiv must be equal to the number of inputs');
                end
                if ~isfield(ser, 'rnamesi')
                    error('There must be option rnamesi with option rnamesiv')
                elseif ser.rnamesi ~= 1
                    error('Value of rnamesi must be equal to one if rnamesiv is also option');
                end
                ser.rnamesiv = val{i};
            end
        elseif strcmpi('prelivar', opt{i})
            if (~isscalar(val{i}) || (val{i} ~= 0 && val{i} ~= 1))
                error('Value of option prelivar must be zero or one');
            else
                ser.prelivar = val{i};
            end
        elseif strcmpi('tfident', opt{i})
            if (~isscalar(val{i}) || (val{i} ~= 0 && val{i} ~= 1))
                error('Value of option tfident must be zero or one');
            else
                ser.tfident = val{i};
            end
        elseif strcmpi('delay', opt{i})
            if ~isnumeric(val{i})
                error('Value of option delay must be a numeric (1 x ninput) array');
            else
                if size(val{i}, 2) ~= size(x, 2)
                    error('Number of delay must be equal to the number of inputs');
                end
                ser.delay = val{i};
            end
        elseif strcmpi('ma', opt{i})
            if ~isnumeric(val{i})
                error('Value of option ma must be a numeric (1 x ninput) array');
            else
                if size(val{i}, 2) ~= size(x, 2)
                    error('Number of ma must be equal to the number of inputs');
                end
                ser.ma = val{i};
            end
        elseif strcmpi('ar', opt{i})
            if ~isnumeric(val{i})
                error('Value of option ar must be a numeric (1 x ninput) array');
            else
                if size(val{i}, 2) ~= size(x, 2)
                    error('Number of ar must be equal to the number of inputs');
                end
                ser.ar = val{i};
            end
        elseif strcmpi('inc', opt{i})
            if (~isscalar(val{i}) || (val{i} ~= 0 && val{i} ~= 1))
                error('Value of option inc must be zero or one');
            else
                ser.inc = val{i};
            end
        elseif strcmpi('modinput', opt{i})
            if ~isstruct(val{i})
                error('Value of option modinput must be a (1 x ninput) structure');
            else
                ser.modinput = val{i};
            end
        elseif strcmpi('modpred', opt{i})
            if ~isstruct(val{i})
                error('Value of option modpred must be a (1 x ninput) structure');
            else
                ser.modpred = val{i};
            end
        elseif strcmpi('backwd', opt{i})
            if (~isscalar(val{i}) || (val{i} ~= 0 && val{i} ~= 1))
                error('Value of option backwd must be zero or one');
            else
                ser.backwd = val{i};
            end
        elseif strcmpi('Cb', opt{i})
            if ~isscalar(val{i})
                error('Value of option Cb must be a scalar');
            else
                ser.Cb = val{i};
            end
        elseif strcmpi('nlagtf', opt{i})
            if (~isscalar(val{i}) || (val{i} < 0 && val{i} ~= -1))
                error('Value of option nlagtf must be positive or minus one');
            else
                ser.nlagtf = val{i};
            end
        elseif strcmpi('maxndtf', opt{i})
            if (~isscalar(val{i}) || (val{i} < 0))
                error('Value of option maxndtf must be positive');
            else
                ser.maxndtf = val{i};
            end
        elseif strcmpi('maxddtf', opt{i})
            if (~isscalar(val{i}) || (val{i} < 0))
                error('Value of option maxddtf must be positive');
            else
                ser.maxddtf = val{i};
            end
        end
    end
end

%call arimaestwi to handle the tf model according to the instructions
%provided by the user
%printing will be done in the subdirectory results. We create it if it does
%not exits
pathc = pwd; %current directory
if exist([pathc, filesep, 'results'], 'dir') ~= 7
    mkdir(pathc, 'results'); %create directory results
end
if isfield(ser, 'gft')
    gft = ser.gft;
    if (gft > 1)
        if exist([pathc, filesep, 'graphs'], 'dir') ~= 7
            mkdir(pathc, 'graphs'); %create directory graphs
        end
    end
else
    gft = 0;
end

out = arimaestwi(sname, ser);

return

%save figures
if (gft > 1)
    datei = cal(ser.bg_year, ser.bg_per, freq);
    if isfield(out.model, 'nmiss')
        nmiss = out.model.nmiss;
        yinterp = out.model.yinterp;
    else
        nmiss = 0;
    end
    s = out.freq;
    lam = out.model.lam;
    %replace missing values, if any, with tentative values as in arimaestni
    y = chmarima(y);
    if lam == 0
        yor = y;
        y = log(yor);
    else
        yor = y;
    end
    ny = length(yor);
    cw = 1.96;
    if isfield(out.model, 'npr')
        npr = out.model.npr;
        pry = out.model.pry;
        spry = out.model.spry;
        if lam == 0
            opry = out.model.opry;
        end
    else
        npr = 0;
    end
    dbname = strrep(sname, '_', '\_');
    %Youtg: outlier effects, Yrg: regression effects.
    %   plotres(y,Y,g,yor,datei,cw,dbname,gf11,nrout,Youtg,mreg,Yrg,infr,s,lam);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    if isfield(out.model, 'Y')
        Y = out.model.Y;
        g = out.model.hb;
    else
        Y = [];
        g = [];
    end
    infr = out.tfmodel.resinf;
    fname = dbname;
    Youtg = out.model.Youtg;
    if isempty(Youtg)
        nrout = 0;
    else
        nrout = size(Youtg, 2);
    end
    Yrg = out.model.Yrg;
    if isempty(Yrg)
        nreg = 0;
    else
        nreg = size(Yrg, 2);
    end
    if gft == 3
        gflag = 1;
    else
        gflag = 0;
    end
    if isfield(infr, 'ycii')
        ycii = infr.ycii;
    else
        ycii = [];
    end
    
    e = infr.e;
    ne = infr.ne;
    stde = infr.stde;
    r = infr.r;
    pc = infr.pc;
    sea = infr.sea;
    sep = infr.sep;
    no = infr.no;
    xo = infr.xo;
    rs = infr.rs;
    pcs = infr.pcs;
    seas = infr.seas;
    
    if ~isempty(y)
        ny = length(y);
        if gft == 2
            f = figure('visible', 'off');
        else
            f = figure('visible', 'on');
        end
        if (s == 12) || (s == 4)
            tsplot(yor, datei, fname)
        else
            plot(yor);
            legend(fname);
        end
        saveas(f, [pathc, filesep, 'graphs', filesep, 'Original'], 'pdf')
        if gflag == 1, disp('strike any key when ready');
            pause;
        end
        flagycii = 0;
        if ~isempty(ycii)
            if gft == 2
                f = figure('visible', 'off');
            else
                f = figure('visible', 'on');
            end
            if lam == 0
                vnames = char('Original Series in logs','Series corrected by filtered inputs only');
            else
                vnames = char('Original Series','Series corrected by filtered inputs only');
            end
            if (s == 12) || (s == 4)
                tsplot([y, ycii], datei, vnames);
            else
                t = 1:ny;
                plot(t, y, t, ycii);
                legend(vnames);
            end
            saveas(f, [pathc, filesep, 'graphs', filesep, 'Originalcorrfi'], 'pdf')
            if gflag == 1, disp('strike any key when ready');
                pause;
            end
            flagycii = 1;
        end
        [ng, mg] = size(g);
        if ng > 0
            if gft == 2
                f = figure('visible', 'off');
            else
                f = figure('visible', 'on');
            end
            if flagycii == 0
                ycor = y - Y(1:ny, :) * g;
            else
                ycor = ycii - Y(1:ny, :) * g;
            end
            if lam == 0
                vnames = char('Original Series in logs','Corrected Series');
            else
                vnames = char('Original Series','Corrected Series');
            end
            if (s == 12) || (s == 4)
                tsplot([y, ycor], datei, vnames);
            else
                t = 1:ny;
                plot(t, y, t, ycor);
                legend(vnames);
            end
            saveas(f, [pathc, filesep, 'graphs', filesep, 'Originalcorr'], 'pdf')
            if gflag == 1, disp('strike any key when ready');
                pause;
            end
        end
    end
    
    if nreg > 0
        %produce graph of regression effects other than outliers
        if gft == 2
            f = figure('visible', 'off');
        else
            f = figure('visible', 'on');
        end
        [nYr, mYr] = size(Yrg);
        t = 1:nYr;
        Yrgt = zeros(nYr, 1);
        for i = 1:nYr
            Yrgt(i) = sum(Yrg(i, 1:nreg));
        end
        plot(t, Yrgt);
        title('Sum of regression effects other than outliers'); %hold off;
        saveas(f, [pathc, filesep, 'graphs', filesep, 'Sumofregef'], 'pdf')
        if gflag == 1
            disp('strike any key when ready')
            pause
        end
    end
    
    if nrout > 0
        %produce graph of outlier effects
        if gft == 2
            f = figure('visible', 'off');
        else
            f = figure('visible', 'on');
        end
        [nYo, mYo] = size(Youtg);
        t = 1:nYo;
        for i = 1:nrout
            plot(t, Youtg(:, i));
            hold on
        end
        title('Outlier effects');
        hold off;
        saveas(f, [pathc, filesep, 'graphs', filesep, 'Outliereff'], 'pdf')
        if gflag == 1
            disp('strike any key when ready')
            pause
        end
    end
    
    if gft == 2
        f = figure('visible', 'off');
    else
        f = figure('visible', 'on');
    end
    t = 1:ne;
    bd = ones(ne, 1) * (cw * stde); %confidence bands
    plot(t, e, t, zeros(ne, 1), t, bd, 'r', t, -bd, 'r') %plot residuals
    legend('Residuals')
    saveas(f, [pathc, filesep, 'graphs', filesep, 'Residuals'], 'pdf')
    if gflag == 1
        disp('strike any key when ready')
        pause
    end
    
    if gft == 2
        f = figure('visible', 'off');
    else
        f = figure('visible', 'on');
    end
    bar(xo, no);
    title('Residual histogram')
    xlabel('Standard deviation intervals')
    saveas(f, [pathc, filesep, 'graphs', filesep, 'Residualhist'], 'pdf')
    if gflag == 1
        disp('strike any key when ready')
        pause
    end
    
    if gft == 2
        f = figure('visible', 'off');
    else
        f = figure('visible', 'on');
    end
    fname = 'residuals';
    rpplot(r, pc, sea, sep, cw, fname); %plot residual autocorrelations
    saveas(f, [pathc, filesep, 'graphs', filesep, 'Resoutcor'], 'pdf')
    if gflag == 1
        disp('strike any key when ready')
        pause
    end
    
    if gft == 2
        f = figure('visible', 'off');
    else
        f = figure('visible', 'on');
    end
    fname = 'squared residuals';
    rpplot(rs, pcs, seas, sep, cw, fname); %plot autocorrelations
    saveas(f, [pathc, filesep, 'graphs', filesep, 'SqResoutcor'], 'pdf')
    if gflag == 1
        disp('strike any key when ready')
        pause
    end
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if (nmiss > 0)
        if gft == 2
            f = figure('visible', 'off');
        else
            f = figure('visible', 'on');
        end
        vnames = char('Interpolated series');
        tsplot(yinterp, datei, vnames);
        saveas(f, [pathc, filesep, 'graphs', filesep, 'interp'], 'pdf')
        if gft == 3
            disp('press any key to continue');
            pause;
        end
    end
    if npr > 0
        %plot forecasts
        if gft == 2
            f = figure('visible', 'off');
        else
            f = figure('visible', 'on');
        end
        tt = ny - npr + 1:ny;
        y1 = [y(tt); pry];
        y2 = [y(tt); pry + cw * spry];
        y3 = [y(tt); pry - cw * spry];
        t = -npr + 1:npr;
        vnames = char('Upper 95% band','Forecast','Lower 95% band');
        plot(t, y2, '-.', t, y1, 'r-', t, y3, '--');
        legend(vnames);
        set(gca, 'tickdir', 'in');
        set(gca, 'xcolor', 'k');
        set(gca, 'GridLineStyle', ':');
        set(gca, 'Xgrid', 'on');
        axis tight;
        if lam == 0
            title(['Forecasts of series ', dbname, ' (in logs)'])
            saveas(f, [pathc, filesep, 'graphs', filesep, 'forecastl'], 'pdf')
            if gft == 3
                disp('press any key to continue');
                pause;
            end
            if gft == 2
                f = figure('visible', 'off');
            else
                f = figure('visible', 'on');
            end
            y1 = [yor(tt); opry];
            t = -npr + 1:npr;
            plot(t, y1, 'r-');
            set(gca, 'tickdir', 'in');
            set(gca, 'xcolor', 'k');
            set(gca, 'GridLineStyle', ':');
            set(gca, 'Xgrid', 'on');
            axis tight;
            title(['Forecasts of series ', dbname])
            saveas(f, [pathc, filesep, 'graphs', filesep, 'forecast'], 'pdf')
        else
            title(['Forecasts of series  ', dbname])
            saveas(f, [pathc, filesep, 'graphs', filesep, 'forecast'], 'pdf')
        end
        if gft == 3
            disp('press any key to continue');
            pause;
        end
        if (s > 1)
            %plot forecasts in rates of grothw
            if gft == 2
                f = figure('visible', 'off');
            else
                f = figure('visible', 'on');
            end
            yt = tasa([yor; opry], s) * 100;
            nyt = length(yt);
            tt = nyt - 2 * npr + 1:nyt;
            t = -npr + 1:npr;
            plot(t, yt(tt), 'r-');
            set(gca, 'tickdir', 'in');
            set(gca, 'xcolor', 'k');
            set(gca, 'GridLineStyle', ':');
            set(gca, 'Xgrid', 'on');
            axis tight;
            title(['Forecasts of original series ', dbname, ' (rates of growth in percentage)'])
            saveas(f, [pathc, filesep, 'graphs', filesep, 'forecastrg'], 'pdf')
        end
    end
    if gft == 3
        disp('End of series. Save figures before continuing');
        disp('press any key to continue');
        pause;
    end
    close all
end
