function [X, Z, G, W, T, H, ins, ii, strc, ferror] = sucdm(comp, y, Y, stra, npr)
%************************************************************************
% PURPOSE: this function sets up a state space model given a structure
% containing information about the components of a canonical decompostion
% model. It returns a structure containing the model information. The
% canonical decomposition model is
%
%  y_t = Y_t*beta + p_t + s_t + i_t,
%
% where Y_t is a vector of regression variables, p_t is the canonical
% trend-cycle, s_t is the canonical seassonal and i_t is the canonical
% irregular component.
%---------------------------------------------------
% USAGE: [X,Z,G,W,T,H,ins,ii,strc,ferror] = sucdm(comp,y,Y,stra,npr)
% Inputs:
%           comp   = a structure with the following fields:
%          .ptnum  = a polynomial containing the trend-cycle numerator
%          .ptden  = a polynomial containing the trend-cycle denominator
%          .ptnur  = number of unit roots in ptden
%          .ptvar  = a positive number containing the variance of the
%                    trend-cycle innovations
%          .stnum  = a polynomial containing the seasonal numerator
%          .stden  = a polynomial containing the seasonal denominator
%          .stnur  = number of unit roots in stden
%          .stvar  = a positive number containing the variance of the
%                    seasonal innovations
%          .rt     = a polynomial containing the transitory component
%          .rtvar  = a positive number containing the variance of the
%                    transitory component innovations.
%          .itvar  = a positive number containing the variance of the
%                    irregular component
%          .sigmaa = a positive number containing the variance of the
%                    series model innovations.
%           y      = an (n x 1) matrix containing the data
%           Y      = an (n x nbeta) matrix containing the regression
%                    variables
%           stra   = a structure, given as output by function suvarmapqPQ
%           npr    = number of forecasts
%
% Outputs:
%       X,Z,G,W,T,H,ins,ii are the matrices and initial conditions
%       information corresponding to the state space model
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
%       More specifically,
%       X    : a  (p x nbeta) matrix
%              it can be []
%       W    : an empty matrix
%       Z    : a  (p x nalpha) matrix
%       G    : a  (p x nepsilon) matrix
%       T    : an (nalpha x nalpha) matrix
%       H    : an (nalpha x nepsilon) matrix
%       ins: an nalpha x (cc+cw0+ca1+cca1) matrix containing the initial
%            state information, according to array i below
%       ii   : a  1 x 4 array containing 4 integers, ii=[cc cw0 ca1 cca1],
%            where
%            cc   = nalpha if c is not missing (0 if c missing)
%            cw0  = number of columns in W_0 (0 if W_0 missing)
%            ca1  = 1 if a_1 is not missing (0 if a_1 missing)
%            cca1 = number of columns in A_1 (0 if A_1 missing)
%
% Other outputs:
%       strc     = a structure with the following fields:
%       .stra
%       .X,.Z,.G,.W,.T,.H,.ins and .ii
%       .comp
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
    disp('field ptnum must be present in structure comp in sucdm.m')
    ferror = 1;
    return
else
    ptnum = comp.ptnum;
end
if ~isfield(comp, 'ptden')
    disp('field ptden must be present in structure comp in sucdm.m')
    ferror = 1;
    return
else
    ptden = comp.ptden;
end
if ~isfield(comp, 'ptnur')
    disp('field ptnur must be present in structure comp in sucdm.m')
    ferror = 1;
    return
else
    ptnur = comp.ptnur;
end
if ~isfield(comp, 'ptvar')
    disp('field ptvar must be present in structure comp in sucdm.m')
    ferror = 1;
    return
else
    ptvar = comp.ptvar;
end
if ~isfield(comp, 'stnum')
    disp('field stnum must be present in structure comp in sucdm.m')
    ferror = 1;
    return
else
    stnum = comp.stnum;
end
if ~isfield(comp, 'stden')
    disp('field stden must be present in structure comp in sucdm.m')
    ferror = 1;
    return
else
    stden = comp.stden;
end
if ~isfield(comp, 'stnur')
    disp('field stnur must be present in structure comp in sucdm.m')
    ferror = 1;
    return
else
    stnur = comp.stnur;
end
if ~isfield(comp, 'stvar')
    disp('field stvar must be present in structure comp in sucdm.m')
    ferror = 1;
    return
else
    stvar = comp.stvar;
end
if ~isfield(comp, 'rt')
    disp('field rt must be present in structure comp in sucdm.m')
    ferror = 1;
    return
else
    rt = comp.rt;
end
if ~isfield(comp, 'rtvar')
    disp('field rtvar must be present in structure comp in sucdm.m')
    ferror = 1;
    return
else
    rtvar = comp.rtvar;
end
if ~isfield(comp, 'itvar')
    disp('field itvar must be present in structure comp in sucdm.m')
    ferror = 1;
    return
else
    itvar = comp.itvar;
end
if ~isfield(comp, 'sigmaa')
    disp('field sigmaa must be present in structure comp in sucdm.m')
    ferror = 1;
    return
else
    sigmaa = comp.sigmaa;
end


np = length(ptnum);
dp = length(ptden);
if np ~= dp
    disp('ptnum and ptden should have the same length in sucdm')
    ferror = 2;
    return
end
ns = length(stnum);
ds = length(stden);
if ns ~= ds
    disp('stnum and stden should have the same length in sucdm')
    ferror = 2;
    return
end


nr = length(rt);
nalpha = np + ns + nr;
npm1 = np - 1;
nsm1 = ns - 1;
nrm1 = nr - 1;
nppns = np + ns;
% matrix T
Tp = zeros(np);
Tp(1:npm1, 2:end) = eye(npm1);
Tp(np, 2:np) = -ptden(1:end-1);
Ts = zeros(ns);
Ts(1:nsm1, 2:end) = eye(nsm1);
Ts(ns, 2:ns) = -stden(1:end-1);
Tr = zeros(nr);
if nr > 1
    Tr(1:nrm1, 2:end) = eye(nrm1);
end
T = zeros(nalpha);
T(1:np, 1:np) = Tp;
T(np+1:nppns, np+1:nppns) = Ts;
if nr > 0
    nppnspnr = nppns + nr;
    T(nppns+1:nppnspnr, nppns+1:nppnspnr) = Tr;
end
% matrix H
if (nr > 1)
    nch = 4;
else
    nch = 3;
end
H = zeros(nalpha, nch);
psip = poldiv(fliplr(ptnum), fliplr(ptden), npm1);
psis = poldiv(fliplr(stnum), fliplr(stden), nsm1);
psip = psip * sqrt(ptvar) * (sigmaa);
psis = psis * sqrt(stvar) * (sigmaa);
H(1:np, 1) = psip';
H(np+1:nppns, 2) = psis';
if (nr > 0)
    psir = fliplr(rt);
    psir = psir * sqrt(rtvar) * (sigmaa);
    H(nppns+1:end, 3) = psir';
end
% matrix G
G = zeros(1, nch);
G(nch) = sqrt(itvar) * (sigmaa);
% matrices Z, W and X
Z = zeros(1, nalpha);
Z(1) = 1.;
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
ndelta = ptnur + stnur;
[ins, ii, ferror] = incossm(T, H, ndelta);

strc.stra = stra;

strc.X = X;
strc.Z = Z;
strc.G = G;
strc.W = W;
strc.T = T;
strc.H = H;
strc.ins = ins;
strc.i = ii;

strc.compcd = comp;
