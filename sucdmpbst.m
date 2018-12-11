function [X, Z, G, W, T, H, ins, ii, strc, ferror] = sucdmpbst(comp, compf, y, Y, ...
    stra, npr)
%************************************************************************
% PURPOSE: this function sets up a state space model given a structure
% containing information about the components of a canonical decompostion
% model, where the trend-cycle is further decomposed into a smooth trend
% and a cycle by means of the application of a low-pass filter of the
% Butterworth sine or tangent type. See "The Use of Butterworth Filters for
% Trend and Cycle Estimation in Economic Time Series", Gómez, V. (2001),
% Journal of Business and Economic Statistics, 19, 365-373. It returns a
% structure containing the model information. The canonical decomposition
% model is
%
%  y_t = Y_t*beta + p_t + s_t + i_t,
%
% where Y_t is a vector of regression variables, p_t is the canonical
% trend-cycle, s_t is the canonical seassonal and i_t is the canonical
% irregular component. The trend-cycle is further decomposed as
%
% p_t = ps_t + c_t,
%
% where the models of ps_t and c_t are given by the parameters of p_t and
% the parameters of the filter.
%---------------------------------------------------
% USAGE: [X,Z,G,W,T,H,ins,ii,strc,ferror] = sucdmpbst(comp,compf,y,Y,stra,npr)
% Inputs:
%             comp  = a structure with the following fields:
%           .ptnum  = a polynomial containing the trend-cycle numerator
%           .ptden  = a polynomial containing the trend-cycle denominator
%           .ptnur  = number of unit roots in ptden
%           .ptvar  = a positive number containing the variance of the
%                     trend-cycle innovations
%           .stnum  = a polynomial containing the seasonal numerator
%           .stden  = a polynomial containing the seasonal denominator
%           .stnur  = number of unit roots in stden
%           .stvar  = a positive number containing the variance of the
%                     seasonal innovations
%           .rt     = a polynomial containing the transitory component
%           .rtvar  = a positive number containing the variance of the
%                    transitory component innovations.
%           .itvar  = a positive number containing the variance of the
%                     irregular component
%           .sigmaa = a positive number containing the variance of the
%                     series model innovations.
%           .phi    = a polynomial containing the regular phi
%
%          compf    = a structure with the following fields
%          .num     = a polynomial containing the filter numerator
%          .den     = a polynomial containing the filter denominator
%          .Alpha   = a polynomial containing the filter Alpha
%          .sa      = a number, filter sa
%          .Di      = a positive integer, filter Di
%          .Thetac  = a number, the frequency, divided by pi, of gain .5 in
%                     the But. sine/tangent filter.
%          .Lambda  = a positive number, filter Lambda
%
%           y       = an (n x 1) matrix containing the data
%           Y       = an (n x nbeta) matrix containing the regression
%                     variables
%           stra    = a structure, given as output by function suvarmapqPQ
%           npr     = number of forecasts
%
% Outputs:
%      X,Z,G,W,T,H,ins,ii are the matrices and initial conditions
%      information corresponding to the state space model
%
%       y_t = X_t*beta + Z_t*alpha_t + G_t*epsilon_t
%       alpha_{t+1}= W_t*beta + T_t*alpha_t + H_t*epsilon_t,
%
%       where epsilon_t is (0,sigma^2I),
%
%       with initial state
%
%       alpha_1= c + W_0*beta + a_1 + A_1*delta
%
%       where c is (0,Omega) and delta is (0,kI) (diffuse).
%
% More specifically,
%       X    : a  (p x nbeta) matrix
%              it can be []
%       W    : an empty matrix
%       Z    : a  (p x nalpha) matrix
%       G    : a  (p x nepsilon) matrix
%       T    : an (nalpha x nalpha) matrix
%       H    : an (nalpha x nepsilon) matrix
%       ins: an nalpha x (cc+cw0+ca1+cca1) matrix containing the initial
%            state information, according to array i below
%       ii   : a  1 x 4 array containing 4 integers, i=[cc cw0 ca1 cca1],
%            where
%            cc   = nalpha if c is not missing (0 if c missing)
%            cw0  = number of columns in W_0 (0 if W_0 missing)
%            ca1  = 1 if a_1 is not missing (0 if a_1 missing)
%            cca1 = number of columns in A_1 (0 if A_1 missing)
%
% Other outputs:
%        strc     = a structure with the following fields:
%        .stra
%        .X,.Z,.G,.W,.T,.H,.ins and .ii
%        .comp
%        .compf
%
%
% Copyright (c) 21 July 2003 by Victor Gomez
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

if ~isfield(comp, 'ptnum')
    disp('field ptnum must be present in structure comp in sucdmpbst.m')
    ferror = 1;
    return
else
    ptnum = comp.ptnum;
end

if ~isfield(comp, 'ptden')
    disp('field ptden must be present in structure comp in sucdmpbst.m')
    ferror = 1;
    return
else
    ptden = comp.ptden;
end
if ~isfield(comp, 'ptnur')
    disp('field ptnur must be present in structure comp in sucdmpbst.m')
    ferror = 1;
    return
else
    ptnur = comp.ptnur;
end
if ~isfield(comp, 'ptvar')
    disp('field ptvar must be present in structure comp in sucdmpbst.m')
    ferror = 1;
    return
else
    ptvar = comp.ptvar;
end
if ~isfield(comp, 'stnum')
    disp('field stnum must be present in structure comp in sucdmpbst.m')
    ferror = 1;
    return
else
    stnum = comp.stnum;
end
if ~isfield(comp, 'stden')
    disp('field stden must be present in structure comp in sucdmpbst.m')
    ferror = 1;
    return
else
    stden = comp.stden;
end
if ~isfield(comp, 'stnur')
    disp('field stnur must be present in structure comp in sucdmpbst.m')
    ferror = 1;
    return
else
    stnur = comp.stnur;
end
if ~isfield(comp, 'stvar')
    disp('field stvar must be present in structure comp in sucdmpbst.m')
    ferror = 1;
    return
else
    stvar = comp.stvar;
end
if ~isfield(comp, 'rt')
    disp('field rt must be present in structure comp in sucdmpbst.m')
    ferror = 1;
    return
else
    rt = comp.rt;
end
if ~isfield(comp, 'rtvar')
    disp('field rtvar must be present in structure comp in sucdmpbst.m')
    ferror = 1;
    return
else
    rtvar = comp.rtvar;
end
if ~isfield(comp, 'itvar')
    disp('field itvar must be present in structure comp in sucdmpbst.m')
    ferror = 1;
    return
else
    itvar = comp.itvar;
end
if ~isfield(comp, 'sigmaa')
    disp('field sigmaa must be present in structure comp in sucdmpbst.m')
    ferror = 1;
    return
else
    sigmaa = comp.sigmaa;
end
if ~isfield(comp, 'phi')
    disp('field phi must be present in structure comp in sucdmpbst.m')
    ferror = 1;
    return
else
    phi = comp.phi;
end
if ~isfield(compf, 'num')
    disp('field num must be present in structure comp in sucdmpbst.m')
    ferror = 1;
    return
else
    num = compf.num;
end
if ~isfield(compf, 'den')
    disp('field den must be present in structure comp in sucdmpbst.m')
    ferror = 1;
    return
else
    den = compf.den;
end
if ~isfield(compf, 'Di')
    disp('field Di must be present in structure comp in sucdmpbst.m')
    ferror = 1;
    return
else
    Di = compf.Di;
end
if ~isfield(compf, 'Alpha')
    disp('field Alpha must be present in structure comp in sucdmpbst.m')
    ferror = 1;
    return
else
    Alpha = compf.Alpha;
end
if ~isfield(compf, 'sa')
    disp('field sa must be present in structure comp in sucdmpbst.m')
    ferror = 1;
    return
else
    sa = compf.sa;
end
% if ~isfield(compf,'Thetac')
%     disp('field Thetac must be present in structure comp in sucdmpbst.m')
%     ferror=1;
%     return
% else
%  Thetac=comp.Thetac;
% end
if ~isfield(compf, 'Lambda')
    disp('field Lambda must be present in structure comp in sucdmpbst.m')
    ferror = 1;
    return
else
    Lambda = compf.Lambda;
end


np = length(ptnum);
dp = length(ptden);
if np ~= dp
    disp('ptnum and ptden should have the same length in sucdmpbst.m')
    ferror = 2;
    return
end
ns = length(stnum);
ds = length(stden);
if ns ~= ds
    disp('stnum and stden should have the same length in sucdmpbst.m')
    ferror = 2;
    return
end

ndelta = ptnur + stnur;

%
% Cascade implementation in sp_t and c_t for numerical accuracy. See
% function cascade
%
% model for sp_t: phisp(z)sp_t = thsp(z) a_spt, where phisy(z) =
% ptden(z) den(z) and thst(z) = ptnum(z)num(z). Implemented
% as the cascade [num(z)/den(z)][ptnum(z)/ptden(z)]. There is no
% cancellation of terms.
nden = length(den);
nnum = length(num);
if nnum < nden
    numx = zeros(size(den));
    numx(end) = 1.; %Butterworth sine
else
    numx = num; %Butterworth tangent
end
sigmap = sqrt(ptvar) * (sigmaa) * (1 ./ sa);
[Tsp, Hsp, Zsp, ferror] = cascade(den, numx, ptden, ptnum, sigmap);
% model for c_t: phic(z) c_t = thc(z) a_ct, where phic(z) =
% phi(z)den(z) and thc(z) = ptnum(z)(1-z)^(Di-ptnur). Implemented
% if Di <= ptnur as the cascade [Alphac(z)/denc(z)][1.], where
% Alphac(z)=ptnum(z) and denc=den*(1-z)^(ptnur-Di). If Di > ptnur,
% implemented as the cascade [Alphac(z)/den(z)][(1-z)^p/phi(z)], where
% Alphac(z)=ptnum(z)(1-z)^(Di-ptnur-p) and p denotes the degree of
% the polynomial phi(z). The factor (1-z)^(ptnur - Di) is
% cancelled.
Alphac = Alpha;
mptndi = min(ptnur, Di);
for i = 1:mptndi
    Alphac = deconv(Alphac, [-1., 1.]);
end
if Di < ptnur
    disp('The number of unit roots in the ARIMA model should be')
    disp('equal or less than Di in the filter. Otherwise, there')
    disp('will be numerical problems')
    ferror = 1;
    return
elseif Di == ptnur
    den = conv(den, phi);
    phi = [0., 1.];
    S = [0., 1];
else
    p = length(phi) - 1;
    S = 1.;
    if p == 0
        phi = [0., 1.];
        S = [0., 1];
    else
        for i = 1:p
            Alphac = deconv(Alphac, [-1., 1.]);
            S = conv(S, [-1., 1.]);
        end
    end
end
Alphac = conv(Alphac, ptnum);

sigmac = sqrt(ptvar) * (sigmaa) * (Lambda / sa);
[Tc, Hc, Zc, ferror] = cascade(den, Alphac, phi, S, sigmac);

% matrix T
[nsp, junk] = size(Tsp);
[nc, junk] = size(Tc);
np = nsp + nc;
nr = length(rt);
nrm1 = nr - 1;
nalpha = np + ns + nr;
nppns = np + ns;
T = zeros(nalpha);
T(1:nsp, 1:nsp) = Tsp;
T(nsp+1:nsp+nc, nsp+1:nsp+nc) = Tc;
Ts = zeros(ns);
nsm1 = ns - 1;
Ts(1:nsm1, 2:end) = eye(nsm1);
Ts(ns, 2:ns) = -stden(1:end-1);
Tr = zeros(nr);
if nr > 1
    Tr(1:nrm1, 2:end) = eye(nrm1);
end
T(np+1:nppns, np+1:nppns) = Ts;
if nr > 0
    nppnspnr = nppns + nr;
    T(nppns+1:nppnspnr, nppns+1:nppnspnr) = Tr;
end
% matrix H
if (nr > 1)
    nch = 5;
else
    nch = 4;
end
H = zeros(nalpha, nch);
H(1:nsp, 1) = Hsp;
H(nsp+1:nsp+nc, 2) = Hc;
psis = poldiv(fliplr(stnum), fliplr(stden), nsm1);
psis = psis * sqrt(stvar) * (sigmaa);
H(np+1:np+ns, 3) = psis';
if (nr > 0)
    psir = fliplr(rt);
    psir = psir * sqrt(rtvar) * (sigmaa);
    H(nppns+1:end, 4) = psir';
end
% matrix G
G = zeros(1, nch);
G(nch) = sqrt(itvar) * (sigmaa);
% matrices Z, W and X
Z = zeros(1, nalpha);
Z(1) = 1.;
Z(nsp+1) = 1;
Z(np+1) = 1.;
if (nr > 0)
    Z(nppns+1) = 1.;
end
W = [];
ny = size(y, 1);
if isempty(Y)
    X = Y;
else
    X = Y(1:ny+npr, :);
end

% initial conditions: matrices ins and ii
[ins, ii, ferror] = incossm(T, H, ndelta);

strc.X = X;
strc.Z = Z;
strc.G = G;
strc.W = W;
strc.T = T;
strc.H = H;
strc.ins = ins;
strc.i = ii;

strc.stra = stra;
strc.compcd = comp;
strc.compf = compf;

end
