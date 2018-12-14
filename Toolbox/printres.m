function printres(fid, infr)
%*************************************************************************
% This function prints test results based on residuals
%
%   INPUTS:
%     fid : file identifier, needed for writing
%    infr : residuals structure (output of rescomp)
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
%*************************************************************************

if ~isstruct(infr)
    error('printres: requires a residual structure');
end;

e = infr.e;
ne = infr.ne;
ve = infr.ve;
stde = infr.stde;
sconp = infr.sconp;
conp = infr.conp;
orders = infr.orders;
r = infr.r;
pc = infr.pc;
qstat = infr.qstat;
pval = infr.pval;
df = infr.df;
sea = infr.sea;
sep = infr.sep;
hot0 = infr.hot0;
ho = infr.ho;
me = infr.me;
rstd = infr.rstd;
rtval = infr.rtval;
mde = infr.mde;
skew = infr.skew;
kurt = infr.kurt;
bst = infr.bst;
pnt = infr.pnt;
tsk = infr.tsk;
tkr = infr.tkr;
dw = infr.dw;
ptdw = infr.ptdw;
n0 = infr.n0;
n1 = infr.n1;
nr = infr.nr;
Tval = infr.Tval;
rs = infr.rs;
pcs = infr.pcs;
qstats = infr.qstats;
pvals = infr.pvals;
dfs = infr.dfs;
seas = infr.seas;
h = infr.h;
H = infr.H;
pH = infr.pH;
aic = infr.aic;
bic = infr.bic;

% print information criteria
fprintf(fid, '\n');
fprintf(fid, 'Information criteria:');
clear in
in.rnames = char('  ', 'AIC', 'BIC');
in.fmt = char('%15.4f');
in.fid = fid;
r2 = [aic, bic]';
mprint(r2, in);
fprintf(fid, '\n');

ncol = min(length(e), 12);
datei = cal(1900, 1, ncol); %the year is immaterial, freq. is the number of columns
clear inft
inft.fid = fid;
inft.fh = 0;
inft.wd = 9;
inft.nd = 3;
inft.scale = 0;
fprintf(fid, 'Residuals:\n');
tabla(e, datei, inft) %print residuals
fprintf(fid, '\n');
if ~isempty(ho)
    % print outliers
    clear in
    in.cnames = char('Obs', 'T-value');
    in.fmt = char('%4.f', '%10.4f');
    in.fid = fid;
    fprintf(fid, 'Residuals with |t| > 3.5:\n');
    mprint([hot0', ho'], in);
    fprintf(fid, '\n');
end

% print autocorrelations and partial autocorrelations
clear in
in.cnames = char('Autcor', 'SE', 'Q-stats', 'DF', 'P-values', 'Parcor', 'SE');
rnames = 'Order';
for j = 1:length(orders)
    rnames = char(rnames, num2str(j));
end;
in.rnames = rnames;
in.fmt = char('%10.4f', '%10.4f', '%10.4f', ...
    '%2.f', '%8.2f', '%10.4f', '%10.4f');
in.fid = fid;
fprintf(fid, 'Sample autocorrelations and partial autocorrelations:\n');
mprint([r, sea, qstat, df, pval, pc, sep], in);
fprintf(fid, '_________\n');
fprintf(fid, 'When DF is positive, P-values should be greater than 0.05\n');
fprintf(fid, '(at the 5 per cent level) for model adequacy\n');

%print Pierce Qs value
if isfield(infr, 'Qs')
    Qs = infr.Qs;
    dfQs = infr.dfQs;
    pvalQs = infr.pvalQs;
    nQs = length(Qs);
    fprintf(fid, '\nPierce Qs:\n');
    clear in
    in.cnames = char('Qs-stat', 'DF', 'P-value');
    rnames = 'Order';
    for i = 1:nQs
        rnames = char(rnames, num2str(i));
    end
    in.rnames = rnames;
    in.fmt = char('%10.4f', '%2.f', '%8.2f');
    in.fid = fid;
    TQs = [Qs, dfQs, pvalQs];
    mprint(TQs, in);
    fprintf(fid, '\n');
end

%print Bell's seasonal dummies F-test
if isfield(infr, 'pFtest')
    pF = infr.pFtest;
    F = infr.F;
    dfn = infr.dfn;
    dfd = infr.dfd;
    fprintf(fid, '\n');
    fprintf(fid, 'Stable seasonality F-test based on residuals:\n');
    clear in
    in.cnames = char('DFn', 'DFd', 'stat', 'P-value');
    in.rnames = char('   ', 'F');
    in.fmt = char('%3.f', '%3.f', '%10.4f', '%8.2f');
    in.fid = fid;
    rf = [dfn, dfd, F, pF];
    mprint(rf, in);
    fprintf(fid, '\n');
end


%print residual diagnostics
fprintf(fid, '\n');
fprintf(fid, 'Residual diagnostics:');
clear in
in.rnames = char('  ', 'Sample size ');
in.fmt = char('%10.f');
in.fid = fid;
r0 = [ne];
mprint(r0, in);

clear in
in.rnames = char('  ', 'Median', 'Mean', ...
    'Std of mean', 'T-value of mean');
in.fmt = char('%10.4f');
in.fid = fid;
r1 = [mde, me, rstd, rtval]';
mprint(r1, in);
fprintf(fid, '\n');

clear in
in.cnames = char(' ', 'P-values');
in.rnames = char('   ', 'Normality (BS) ', 'Skewness', 'Kurtosis');
in.fmt = char('%10.4f', '%8.2f');
in.fid = fid;
r2 = [[bst, skew, kurt]', [pnt, tsk, tkr]'];
mprint(r2, in);
fprintf(fid, '\n');

clear in
in.cnames = char(' ', 'P-value');
in.rnames = char('   ', 'Durbin-Watson  ');
in.fmt = char('%10.4f', '%8.2f');
in.fid = fid;
r3 = [dw, ptdw];
mprint(r3, in);

clear in
in.rnames = char('   ', 'Standard error ', 'Sigma square', ...
    'Residual variance', 'Residual std. dev.');
in.fmt = char('%10.4f');
in.fid = fid;
r4 = [sconp, conp, ve, stde]';
mprint(r4, in);
fprintf(fid, '\n');

fprintf(fid, 'Approximate test of runs on residuals:');
clear in
in.rnames = char('  ', 'Number of runs', 'Number of (+)', 'Number of (-)');
in.fmt = char('%8.f');
in.fid = fid;
mprint([nr, n1, n0]', in);

clear in
in.rnames = char('  ', 'T-value      ');
in.fmt = char('%12.4f');
in.fid = fid;
mprint([Tval]', in);

% print heteroscedasticity test
fprintf(fid, '\n');
fprintf(fid, 'Heteroscedasticity test:\n');
clear in
in.cnames = char('DF', 'stat', 'P-value');
in.rnames = char('   ', 'H');
in.fmt = char('%3.f', '%10.4f', '%8.2f');
in.fid = fid;
r2 = [h, H, pH];
mprint(r2, in);


% print autocorrelations and partial autocorrelations of squared residuals
clear in
in.cnames = char('Autcor', 'SE', 'Q-stats', 'DF', 'P-values', 'Parcor', 'SE');
rnames = 'Order';
for j = 1:length(orders);
    rnames = char(rnames, num2str(j));
end;
in.rnames = rnames;
in.fmt = char('%10.4f', '%10.4f', '%10.4f', '%2.f', ...
    '%8.2f', '%10.4f', '%10.4f');
in.fid = fid;
fprintf(fid, '\nSample autocorrelations of squared residuals:\n');
mprint([rs, seas, qstats, dfs, pvals, pcs, sep], in);
fprintf(fid, '_________\n');
fprintf(fid, 'When DF is positive, P-values should be greater than 0.05\n');
fprintf(fid, '(at the 5 per cent level) for model adequacy\n');
%end printing residual diagnostics
