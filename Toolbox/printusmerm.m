function printusmerm(fid, datei, tname, yor, y, ny, lam, modescr, result, nreg, nbeta)
%*************************************************************************
% This function prints the estimation results of a univariate structural
% model with complex seasonal patterns
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
if ~isfield(modescr, 'seas')
    N = 0;
else
    N = length(modescr.seas);
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
%we do not print header and years, only the data. The number of columns is
%given by structure datei in the third field (freq)
% inft = minft(fid,1,9,3,1);
inft = minft(fid, 0, 9, 3, 1);
dateib = datei;
inftb = inft;
fnameb = tname;
fprintf(fid, 'Original series:\n');
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
cc = cell(1, 5+N);
cc(1:3) = {'Sigma.irreg', 'Sigma.level', 'Sigma.slope'};
for i = 1:N
    cc{3+i} = ['Sigma.seasp{', num2str(i), '}'];
end
cc(4+N:5+N) = {'Sigma.autor', 'Sigma.cycle'};
mm = char(cc(conc));
% rnamess=char('Sigma irreg.','Sigma level ','Sigma slope ',...
%          'Sigma seaso.','Sigma autor.','Sigma cycle');
in.rnames = char(in.rnames, mm);
if nst > 0
    for i = 1:nst,
        mm = char(cc(stord(i)));
        in.rnames = char(in.rnames, mm);
    end
end
if cycle > 0
    in.rnames = char(in.rnames, 'Cycle rho ', 'Cycle freq.');
end
if arp > 0
    %  rnamesa=[];
    %  for i=1:arp
    %   rnamesa=[rnamesa; ['Ar(' num2str(i) ')       ']];
    %  end
    for i = 1:arp
        mm = ['Ar(', num2str(i), ')       '];
        in.rnames = char(in.rnames, mm);
    end
end
in.fmt = char('%12.4f', '%12.4f', '%12.4f');
in.fid = fid;
mprint(z, in);
fprintf(fid, ['Parameter ', char(cc(conc)), ' is concentrated out of the likelihood\n']);


Pevf = result.Pevf; % prediction error variance
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
    if flevel == 1, rnames = char(rnames, ['level']);
    end
    if fslope == 1, rnames = char(rnames, ['slope']);
    end
    if nreg > 0
        for i = 1:nreg, rnames = char(rnames, ['reg', num2str(i)]);
        end
    end
    in.rnames = rnames;
    in.fmt = char('%12.5f', '%12.5f', '%8.2f');
    in.fid = fid;
    mprint(Mbeta, in);
    fprintf(fid, '\n');
end
