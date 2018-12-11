function Y = arimasimeasy(freq, varargin)
%**********************************************************************
%                       EASY ARIMA SIMULATION
%
%                            USAGE :
% Y = arimasimeasy(freq,'option1',optionvalue1,'option2',optionvalue2,...)
%
%       INPUTS :
%------------------
%      REQUIRED
%         freq : data frequency (number of observations per year)
%------------------
%       OPTIONS
%              'phir': (1 x p) array containing the regular AR polynomial
%                      in MATLAB format, for example, [-.4 1] for phir = 1.
%                      - .4B, default 1.
%              'phis': (1 x ps) array containing the seasonal AR polynomial
%                      in MATLAB format, default 1.
%               'thr': (1 x q) array containing the regular MA polynomial
%                      in MATLAB format, default 1.
%               'ths': (1 x q) array containing the seasonal MA polynomial
%                      in MATLAB format, default 1.
%          '[p dr q]': (1 x 3) array containing the regular orders
%                      default: [0 0 0]
%        '[ps ds qs]': (1 x 3) array containing the first seasonal orders
%                      default: [0 0 0]
%                 'N': length of the series to be generated, default 100
%                'Ns': number of series to be generated, default 1
%           'discard': number of initial observations to be discarded when
%                      simulating the series, default 50
%               'gft': flag for graphics, = 1, plot series, = 0, no plots
%                       = 2, plots are saved but not displayed, = 3, plots
%                       are both saved and displayed, default 0
%               'drg': regular differencing to be applied to the simulated
%                      series when generating graphs, default 0
%               'dsg': seasonal differencing to be applied to the simulated
%                      series when generating graphs, default 0
%              'mean': mean value for the differenced series, ~=0 mean, =0,
%                       no mean, default 0
%              'stda': standard deviation of the simulated series, default
%                      1.
%              'seed': seed used for the simulation, default 20
%               'lag': number of lags for autocorrelations, default 3*freq
%                'cw': confidence bands coefficient, default 1.96
%
%       OUTPUT : Y   : the simulated (N x Ns) series array
%------------------
%
%     Examples:
%
%   Y=arimasimeasy(freq,'mean',.5)
%   Y=arimasimeasy(freq,'[p dr q]',[0 1 1],'thr',[-.4 1.],'gft',2)
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

if nargin < 1
    error('There must be at least one input to arimaeasy');
end

% set default values
def = cell(1, 18);
ldef = length(def);
for i = 1:4
    def{i} = 1.;
end
def{5} = [0, 0, 0];
def{6} = [0, 0, 0];
def{7} = 100;
def{8} = 1;
def{9} = 50;
def{10} = 1;
def{11} = 0;
def{12} = 0;
def{13} = 0;
def{14} = 1.;
def{15} = 20;
def{16} = max(8, 3*freq);
def{17} = 1.96;
phir = 1.;
phis = 1.;
thr = 1.;
ths = 1.;
p = 0;
dr = 0;
q = 0;
ps = 0;
ds = 0;
qs = 0;
N = 100;
Ns = 1;
discard = 50;
gft = 1;
drg = 0;
dsg = 0;
mean = 0.0;
stda = 1.;
seed = 20;
lag = max(8, 3*freq);
cw = 1.96;

if nargin == 1
    val = def;
end

lvar = length(varargin);

if lvar > 0
    if mod(lvar, 2) == 1 % number of elements is odd
        error('To specify an option, one of the admissible types and an admissible value must be chosen');
    else
        opt = {'phir', ... %1
            'phis', ... %2
            'thr', ... %3
            'ths', ... %4
            '[p dr q]', ... .          %5
            '[ps ds qs]', ... %6
            'N', ... %7
            'Ns', ... %8
            'discard', ... %9
            'gft', ... %10
            'drg', ... %11
            'dsg', ... %12
            'mean', ... %13
            'stda', ... %14
            'seed', ... %15
            'lag', ... %16
            'cw'};%17
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
        if strcmpi('phir', opt{i})
            if (~isnumeric(val{i}) || (size(val{i}, 1) ~= 1))
                error('Value of option phir must be a numeric (1 x p) array');
            else
                phir = val{i};
            end
        elseif strcmpi('phis', opt{i})
            if (~isnumeric(val{i}) || (size(val{i}, 1) ~= 1))
                error('Value of option phis must be a numeric (1 x ps) array');
            else
                phis = val{i};
            end
        elseif strcmpi('thr', opt{i})
            if (~isnumeric(val{i}) || (size(val{i}, 1) ~= 1))
                error('Value of option thr must be a numeric (1 x q) array');
            else
                thr = val{i};
            end
        elseif strcmpi('ths', opt{i})
            if (~isnumeric(val{i}) || (size(val{i}, 1) ~= 1))
                error('Value of option ths must be a numeric (1 x q) array');
            else
                ths = val{i};
            end
        elseif strcmpi('[p dr q]', opt{i})
            if (~isnumeric(val{i}) || (size(val{i}, 1) ~= 1 || size(val{i}, 2) ~= 3))
                error('Value of option [p dr q] must be a numeric (1 x 3) array');
            else
                p = val{i}(1);
                dr = val{i}(2);
                q = val{i}(3);
            end
        elseif strcmpi('[ps ds qs]', opt{i})
            if (~isnumeric(val{i}) || (size(val{i}, 1) ~= 1 || size(val{i}, 2) ~= 3))
                error('Value of option [ps ds qs] must be a numeric (1 x 3) array');
            else
                ps = val{i}(1);
                ds = val{i}(2);
                qs = val{i}(3);
            end
        elseif strcmpi('N', opt{i})
            if (~isscalar(val{i}) || (val{i} <= 0))
                error('Value of option N must be greater than zero');
            else
                N = val{i};
            end
        elseif strcmpi('Ns', opt{i})
            if (~isscalar(val{i}) || (val{i} <= 0))
                error('Value of option Ns must be greater than zero');
            else
                Ns = val{i};
            end
        elseif strcmpi('discard', opt{i})
            if (~isscalar(val{i}) || (val{i} <= 0))
                error('Value of option discard must be greater than zero');
            else
                discard = val{i};
            end
        elseif strcmpi('gft', opt{i})
            if (~isscalar(val{i}) || (val{i} ~= 0 && val{i} ~= 1 && val{i} ~= 2 && val{i} ~= 3))
                error('Value of option gft must be zero, one, two or three');
            else
                gft = val{i};
            end
        elseif strcmpi('drg', opt{i})
            if (~isscalar(val{i}) || (val{i} <= 0))
                error('Value of option drg must be grater than zero');
            else
                drg = val{i};
            end
        elseif strcmpi('dsg', opt{i})
            if (~isscalar(val{i}) || (val{i} <= 0))
                error('Value of option dsg must be grater than zero');
            else
                dsg = val{i};
            end
        elseif strcmpi('mean', opt{i})
            if (~isscalar(val{i}))
                error('Value of option mean must be a scalar');
            else
                mean = val{i};
            end
        elseif strcmpi('stda', opt{i})
            if (~isscalar(val{i}) || (val{i} <= 0))
                error('Value of option stda must be greater than zero');
            else
                stda = val{i};
            end
        elseif strcmpi('seed', opt{i})
            if (~isscalar(val{i}) || (val{i} <= 0))
                error('Value of option seed must be greater than zero');
            else
                seed = val{i};
            end
        elseif strcmpi('lag', opt{i})
            if (~isscalar(val{i}) || (val{i} <= 0))
                error('Value of option lag must be greater than zero');
            else
                lag = val{i};
            end
        elseif strcmpi('cw', opt{i})
            if (~isscalar(val{i}) || (val{i} <= 0))
                error('Value of option cw must be greater than zero');
            else
                cw = val{i};
            end
        end
    end
end

if (p + dr + q + ps + ds + qs) == 0
    disp('model is white noise:')
    disp('simulate it with randn')
    Y = [];
    return
end
if (p > 0) && (length(phir) < p - 1)
    error('length of phir should be p plus 1')
end
if (ps > 0) && (length(phis) < ps - 1)
    error('length of phis should be ps plus 1')
end
if (q > 0) && (length(thr) < q - 1)
    error('length of thr should be p plus 1')
end
if (qs > 0) && (length(ths) < qs - 1)
    error('length of ths should be qs plus 1')
end


%call gensersm to simulate the ARIMA model according to the instructions
%provided by the user
%Create directory graphs if it does not exist when gft > 1
pathc = pwd; %current directory
if (gft > 1)
    if exist([pathc, filesep, 'graphs'], 'dir') ~= 7
        mkdir(pathc, 'graphs'); %create directory graphs
    end
end

% set ARIMA coefficients
x = sarimac(p, ps, q, qs, phir, phis, thr, ths);
%compute ARIMA polynomials
[phi, alpha, th] = arimapol(x, freq, 0, p, ps, dr, ds, 0, q, qs, 0);
% include a mean for the differenced series; cons=0 no mean
if abs(mean) > 0
    ctr = constant(N, mean, dr, ds, 0, 0, freq);
else
    ctr = [];
end

%generate series
Y = gensersm(phi, alpha, th, stda, ctr, Ns, N, discard, seed);

%graphs
if (gft > 0)
    c = acgf(phi, th, lag+1); %theoretical autocovariances
    c0p = c(1);
    cvp = c(2:lag+1);
    rp = cvp / c0p; %theoretical autocorrelations
    [~, pcp] = durlev(c0p, cvp'); %theoretical partial autocorrelations
    if gft == 2
        f = figure('visible', 'off');
    else
        f = figure('visible', 'on');
    end
    arimasplot(Y, dr, drg, freq, dsg, lag, rp, pcp, cw)
    if (gft > 1)
        saveas(f, [pathc, filesep, 'graphs', filesep, 'Simulated'], 'pdf')
    end
    
    if gft ~= 2
        disp('press any key to continue');
        pause;
    end
    close all
end
