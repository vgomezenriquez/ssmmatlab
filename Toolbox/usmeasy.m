function [out, ser] = usmeasy(y, freq, varargin)
%**********************************************************************
%                       EASY STRUCTURAL MODELING
%
%                            USAGE :
% out = usmeasy(y,freq,'option1',optionvalue1,'option2',optionvalue2,...)
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
%                 'Y': (n x nY) array for regression variables, where n is
%                      the series length plus the number of forecasts and
%                      nY is the number of regression variables, default []
%             Ycomp  : a cell array, containing the assignment of each
%                      regression variable to a component. Possible values
%                      are 'level','slope','seas','cycle', 'ar' and 'irreg'
%          'rnamesrg': string matrix for names of regression variables,
%                      default []
%                 'W': (n*nalpha x nbeta) array for the transition equation
%                      of the state space model, where n is the series
%                      length plus the number of forecasts, nalpha is the
%                      state vector length and nbeta is the number of
%                      intervention effects to be modeled this way, default
%                      []
%             'level': (1 x 3) array to specify the level
%             'slope': (1 x 3) array to specify the slope
%             'cycle': (1 x 3) array to specify the cycle
%            'cyclep': (2 x 2) array to specify the rho and alpha
%                      parameters of the cycle
%            'cycleb': (1 x 2) array to specify the cyclical interval
%              'seas': (1 x 3) array to specify the seasonal component
%                'ar': (1 x 3) array to specify the autoregressive
%                      component
%               'arp': (2 x p) array to specify the autoregressive
%                       parameters
%            'conout':'level' if the standard error of the level is
%                      concentrated out
%                      'slope' if the standard error of the slope is
%                      concentrated out
%                      'seas' if the standard error of the seasonal
%                      is concentrated out
%                      'cycle' if the standard error of the cycle is
%                      concentrated out
%                      'ar' if the standard error of the ar component is
%                      concentrated out
%                      'irreg' if the standard error of the irregular is
%                      concentrated out
%                      If .conout is not input, the program will determine
%                      the biggest variance.
%           'sqrtfil': =1 use the square root Kalman filter, =0 do not use
%                      it, default 0
%           'nlestim': flag for nonlinear estimation, = 1, nl estimation,
%                      = 0, no nl estimation, default 1
%               'npr': number of forecasts, default 0
%            'olsres': flag for OLS residuals, = 1, OLS residuals are used,
%                      = 0, uncorrelated residuals (transformation of OLS
%                      residuals) are used, default 0
%                'pr': flag for printing in an external file, = 1, printing
%                       = 0, no printing, default 1
%               'gft': flag for graphics, = 1, plot series, = 0, no plots
%                       = 2, plots are saved but not displayed, = 3, plots
%                       are both saved and displayed, default 0
%             'sname': character array containing the series name
%                      default series1
%-------------------------------------------------------------------------
% * codes for the components:
%   level = -1  constant
%            1  stochastic
%            2  Butterworth tangent
%   slope = -1  constant
%            1  stochastic
%   seas  = -1  fixed dummy seasonality
%            1  stochastic dummy seasonality
%            2  trigonometric seasonality
%            4  Butterworth tangent
%   cycle =  1  structural model cycle
%            2  Butterworth sine cycle
%   irreg =  1  stochastic
%   ar    =  k  autoregressive component of order k
%
%---------------------------------------------------
%
%
%       OUTPUT : a structure, the output of function arimaestni
%------------------
%
%     Example:
%
%   [out,ser]=usmeasy(y,freq,'pr',1,'gft',1,'sname','myusmseries',...
%                     'level',[1 0.1 NaN],'slope',[1 0.1 NaN],'seas',...
%                     [2 0.1 NaN],'irreg',[1 0.1 NaN]);
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
my = size(y, 2);
if my > 1
    error('Array y must be (n x 1) in arimaeasy');
end
if (~isscalar(freq) || (freq < 0))
    error('Value of freq must be a positive integer in arimaeasy');
end


ser.yor = y;
ser.freq = freq;

% set default values for '[bg_year bg_per]'
def = cell(1, 23);
ldef = length(def);
for i = 1:ldef
    def{i} = [];
end
% set default values
def{1} = [2000, 1];
ser.bg_year = 2000;
ser.bg_per = 1;
def{13} = 0;
def{end} = 'usmseries';
sname = 'usmseries';

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
            'Y', ... %3
            'Ycomp', ... %4
            'rnamesrg', ... %5
            'level', ... %6
            'slope', ... %7
            'cycle', ... %8
            'cyclep', ... %9
            'cycleb', ... %10
            'seas', ... %11
            'ar', ... %12
            'arp', ... %13
            'conout', ... %14
            'irreg', ... %15
            'sqrtfil', ... %16
            'nlestim', ... %17
            'npr', ... %18
            'olsres', ... %19
            'pr', ... %20
            'gft', ... %21
            'W', ... %22
            'sname'};%23
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
        elseif strcmpi('Y', opt{i})
            if ~isnumeric(val{i})
                error('Value of option Y must be a numeric array');
            else
                ser.Y = val{i};
            end
        elseif strcmpi('W', opt{i})
            if ~isnumeric(val{i})
                error('Value of option W must be a numeric array');
            else
                ser.W = val{i};
            end
        elseif strcmpi('Ycomp', opt{i})
            if ~iscell(val{i})
                error('Value of option rnamesrg must be a cell array');
            else
                if ~isfield(ser, 'Y')
                    error('There must be option Y with option Ycomp')
                elseif size(val{i}, 2) ~= size(ser.Y, 2)
                    error('Number of Ycomp must be equal to the number of columns of Y');
                end
                ser.Ycomp = val{i};
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
        elseif strcmpi('level', opt{i})
            if (~isnumeric(val{i}) || (size(val{i}, 1) ~= 1 || size(val{i}, 2) ~= 3))
                error('Value of option level must be a numeric (1 x 3) array');
            else
                ser.comp.level = val{i};
            end
        elseif strcmpi('slope', opt{i})
            if (~isnumeric(val{i}) || (size(val{i}, 1) ~= 1 || size(val{i}, 2) ~= 3))
                error('Value of option slope must be a numeric (1 x 3) array');
            else
                ser.comp.slope = val{i};
            end
        elseif strcmpi('cycle', opt{i})
            if (~isnumeric(val{i}) || (size(val{i}, 1) ~= 1 || size(val{i}, 2) ~= 3))
                error('Value of option cycle must be a numeric (1 x 3) array');
            else
                ser.comp.cycle = val{i};
            end
        elseif strcmpi('cyclep', opt{i})
            if ~isfield(ser.comp, 'cycle')
                error('Option cycle should be present when option cyclep is used')
            elseif (~isnumeric(val{i}) || (size(val{i}, 1) ~= 2 || size(val{i}, 2) ~= 2))
                error('Value of option cyclep must be a numeric (2 x 2) array');
            else
                ser.comp.cyclep = val{i};
            end
        elseif strcmpi('cycleb', opt{i})
            if ~isfield(ser.comp, 'cycle') || ~isfield(ser.comp, 'cyclep')
                error('Options cycle and cyclep should be present when option cycleb is used')
            elseif (~isnumeric(val{i}) || (size(val{i}, 1) ~= 1 || size(val{i}, 2) ~= 2))
                error('Value of option cycleb must be a numeric (1 x 2) array');
            else
                ser.comp.cycleb = val{i};
            end
        elseif strcmpi('seas', opt{i})
            if (~isnumeric(val{i}) || (size(val{i}, 1) ~= 1 || size(val{i}, 2) ~= 3))
                error('Value of option seas must be a numeric (1 x 3) array');
            else
                ser.comp.seas = val{i};
            end
        elseif strcmpi('ar', opt{i})
            if (~isnumeric(val{i}) || (size(val{i}, 1) ~= 1 || size(val{i}, 2) ~= 3))
                error('Value of option ar must be a numeric (1 x 3) array');
            else
                ser.comp.ar = val{i};
            end
        elseif strcmpi('arp', opt{i})
            if ~isfield(ser.comp, 'ar')
                error('Option ar should be present when option arp is used')
            elseif (~isnumeric(val{i}) || (size(val{i}, 1) ~= 2 || ...
                    size(val{i}, 2) ~= ser.comp.ar(1)))
                error('Value of option arp must be a numeric (2 x p) array');
            else
                ser.comp.arp = val{i};
            end
        elseif strcmpi('irreg', opt{i})
            if (~isnumeric(val{i}) || (size(val{i}, 1) ~= 1 || size(val{i}, 2) ~= 3))
                error('Value of option irreg must be a numeric (1 x 3) array');
            else
                ser.comp.irreg = val{i};
            end
        elseif strcmpi('conout', opt{i})
            if ~ischar(val{i})
                error('Value of option conout must be a character string');
            else
                ser.comp.conout = val{i};
            end
        elseif strcmpi('sqrtfil', opt{i})
            if (~isscalar(val{i}) || (val{i} ~= 0 && val{i} ~= 1))
                error('Value of option sqrtfil must be zero or one');
            else
                ser.comp.sqrtfil = val{i};
            end
        elseif strcmpi('nlestim', opt{i})
            if (~isscalar(val{i}) || (val{i} ~= 0 && val{i} ~= 1))
                error('Value of option nlestim must be zero or one');
            else
                ser.nlestim = val{i};
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

out = usmestni(sname, ser);
