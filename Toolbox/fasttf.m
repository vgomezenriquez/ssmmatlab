function [f, xv, e, g, M, yd1] = fasttf(xv, y, Y, parm, infm, xf)
%
%
%        this function computes model residuals using the fast Morf, Sihdu
%        and Kailath algorithm
%
%
%        INPUTS:
%       xv: an array containing model parameters
%        y: an array containing the input series
%        Y: a matrix containing regression variables
%     parm: a structure containing model information, where
%       .s:  seasonality
%       .S:  second seasonality
%       .p:  AR order
%      .ps: order of the AR of order s
%       .q:  order of the regular MA
%      .qs: order of the MA of order s (1 at most)
%      .qS: order of the MA of order S (1 at most)
%      .dr: order of regular differencing
%      .ds: order of differencing of order s
%      .dS: order of differencing of order S
%    .pvar:  array containing the indices of variable parameters
%    .pfix:  array containing the indices of fixed parameters
%  .ninput: number of inputs
%  .inputv: array containing the input variables
%   .delay: array with the delays of the input filters
%      .ma: array with the ma parameters of the input filters
%      .ar: array with the ar parameters of the input filters
%     .npr: number of forecasts
%   .nmiss: number of missing values
%     infm: a structure containing information on the estimation method,
%           where
%     .mvx: =1, exact max. likelihood; =0, unconditional mean squares
%     .chb: = 1  compute the beta estimate and its MSE
%             0  do not compute the beta estimate and its MSE
%     .inc: = 0, the initial states in the filter equations to obtain the
%                 filtered variables are equal to zero (not estimated)
%           = 1, the initial states in the filter equations are estimated
%       xf: an array containing fixed parameter values
%
%        OUTPUTS:
%        f: residual vector, whose sum of squares will be minimized
%
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


% parametrization changed 05-2002. ARMA coefficients are used

if ~isstruct(parm)
    error('fasttf: requires the parm structure for optimization');
end;
if ~isstruct(infm)
    error('fasttf: requires the infm structure for optimization');
end;


s = parm.s;
S = parm.S;
dr = parm.dr;
ds = parm.ds;
dS = parm.dS;
p = parm.p;
ps = parm.ps;
q = parm.q;
qs = parm.qs;
qS = parm.qS;
pfix = parm.pfix;
pvar = parm.pvar;
ninput = parm.ninput;
inputv = parm.inputv;
delay = parm.delay;
ma = parm.ma;
ar = parm.ar;
if isfield(parm, 'npr') %npr= number of forecasts
    npr = parm.npr;
else
    npr = 0;
end

mvx = infm.mvx;
chb = infm.chb;
inc = infm.inc;


npar = length(pfix) + length(pvar);
x = zeros(1, npar);
x(pfix) = xf;
x(pvar) = xv;
chk = chkroots(x, p, ps, q, qs, qS);
ppsqqsqS = p + ps + q + qs + qS;
if chk == 1
    xx = invroots(x(1:ppsqqsqS), p, ps, q, qs, qS);
    xx = [xx, x(ppsqqsqS+1:end)];
    xv = xx(pvar);
    x(pvar) = xv;
end
nY = size(Y, 2);
% difference the data
yd = diffest(y, Y, s, S, dr, ds, dS, 0);
nyd = length(yd);

%subtract filtered inputs from output
Yx1 = [];
org = ppsqqsqS;
ny = length(y);
for i = 1:ninput
    inputvd = diffest(inputv(1:ny, i), [], s, S, dr, ds, dS, 0);
    b = delay(i);
    omega = x(org+1:org+ma(i)+1);
    %  if ma(i) > 0
    %   ptf=ma(i); xtf=x(org+2:org+ma(i)+1)/x(org+1);
    %   chk=chkroots(xtf,ptf,0,0,0,0);
    %   if chk == 1
    %    xtf=invroots(xtf,ptf,0,0,0,0);
    %    x(org+2:org+ma(i)+1)=xtf*x(org+1);omega=x(org+1:org+ma(i)+1);
    %    xv=x(pvar);
    %   end
    %   x(pfix)=xf;
    %  end
    delta = 1;
    if ar(i) > 0
        ptf = ar(i);
        omi = org + ma(i);
        xtf = x(omi+2:omi+ptf+1);
        chk = chkroots(xtf, ptf, 0, 0, 0, 0);
        if chk == 1
            xtf = invroots(xtf, ptf, 0, 0, 0, 0);
            x(omi+2:omi+ptf+1) = xtf;
            xv = x(pvar);
        end
        x(pfix) = xf;
        delta = [delta, xtf];
    end
    [z, rx1] = armafil(inputvd, omega, delta, b, inc);
    if isempty(rx1)
        filterisnotconst = 0;
    else
        filterisnotconst = 1;
    end
    %  plot(z)
    %  pause
    org = org + ma(i) + ar(i) + 1;
    yd(:, 1) = yd(:, 1) - z;
    if inc == 1 && filterisnotconst
        %this is the design matrix for the unknown x_1 in the filter
        %equation
        Yx1 = [Yx1, rx1];
    end
end
if (inc == 1) && (ninput > 0) && filterisnotconst
    Yx1 = [Yx1, yd(:, 2:nY+1)];
    %we estimate by regression all of the x_1 and the other regression variables.
    %  beta=Yx1\yd(:,1);
    beta = pinv(Yx1(1:nyd-npr, :)) * yd(1:nyd-npr, 1); %robust OLS estimation
    yd(:, 1) = yd(:, 1) - Yx1(:, 1:end-nY) * beta(1:end-nY);
end

%the following lines added on 3/11/2008
yd1 = yd(:, 1); %store the series corrected by the filtered inputs
if npr > 0
    f = [];
    e = [];
    g = [];
    M = [];
    return
end
%end of lines added on 3/11/2008

% compute ARIMA polynomials for the model
[phi, alprsS, thrsS] = arimapol(x, s, S, p, ps, 0, 0, 0, q, qs, qS);

% set up ARIMA matrices for the model
[~, T, ~, ~, Sigma, ~] = arimam(phi, alprsS, thrsS);

%perform fast likelihood evaluation
if ~isfield(parm, 'nmiss')
    nmiss = 0;
else
    nmiss = parm.nmiss;
end
[F, e, g, M] = fstlkhev(yd1, yd(:, 2:nY+1), T, Sigma, chb, nmiss);
if mvx == 1
    f = F * e;
else
    f = e;
end
