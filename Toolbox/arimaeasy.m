function [out, ser] = arimaeasy(y, freq, varargin)
%**********************************************************************
%                       EASY ARIMA MODELING
%
%                            USAGE :
% out = arimaeasy(y,freq,'option1',optionvalue1,'option2',optionvalue2,...)
%
%       INPUTS :
%------------------
%      REQUIRED
%            y : (ly x 1) array containing the series;
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
%                      default mseries
%
%
%       OUTPUT : a structure, the output of function arimaestni
%------------------
%
%     Examples:
%
%   [out,ser]=arimaeasy(y,freq,'autmid',1,'out',1)
%   out=arimaeasy(y,freq,'[p dr q]',[0 1 1],'leapy',-1)
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

if nargin < 2
    error('There must be at least two inputs to arimaeasy');
end
[ny, my] = size(y);
if my > 1
    error('Array y must be (n x 1) in arimaeasy');
end
if (~isscalar(freq) || (freq < 0))
    error('Value of freq must be a positive integer in arimaeasy');
end


ser.yor = y;
ser.freq = freq;

% set default values for '[bg_year bg_per]'
def = cell(1, 32);
ldef = length(def);
def{1} = [2000, 1];
ser.bg_year = 2000;
ser.bg_per = 1;
for i = 2:ldef - 1
    def{i} = [];
end
% set default value for 'sname'
def{end} = 'mseries';
sname = 'mseries';

if nargin == 2
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
            'sname'};%32
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
        end
    end
end

%call arimaestni to handle the ARIMA model according to the instructions
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
    % else
    %  gft=0;
end

out = arimaestni(sname, ser);
