function infr = rescomp(e, lag, nr, Ss, conp, sconp, Ff, ndrs, nreg)
%*************************************************************************
%        This function generates a structure
%        containing some information on the
%        residuals
%
%    INPUTS:
%        e : residual vector
%      lag : integer specifying lag up to which autocorrelations are to be
%            computed
%       nr : number of parameters to be estimated
%       Ss : residual sum of squares
%     conp : prediction error variance
%    sconp : square root of conp
%       Ff : the product F'*F, where F is the vector of nonlinear functions
%            whose sum of squares is minimized at the end of estimation
%     ndrs : number of residuals for the bias filter e = E*beta + u, where
%            Var(u) = I*sigma^2. It is equal to the number of stacked
%            observations minus the number of nonstationary components and
%            missing values
%     nreg : number of regression variables
%
%   OUTPUT :
%     infr : residual structure containing the following fields:
%       .e : residual vector
%      .ne : length of e
%      .ve : residual variance
%    .stde : residual standard deviation
%    .conp : residual sum of squares
%   .sconp : square root of conp
%  .orders : vector with integers specifying at which lags autocorrelations
%            are to be computed; values between 1 and lag
%       .r : autocorrelations
%      .pc : partial autorrelations
%   .qstat : Q-statistics based on residuals
%    .pval : p-values of the Q-statistics based on residuals
%      .df : degrees of freedom for the Q-statictics based on residuals
%     .sea : standard errors of autocorrelations
%     .sep : standard error of partial correlations;
%            value equal to 1/sqrt(ne)
%      .no : bins of the residual histogram
%      .xo : vector of cut points where observations counted in bin(i)
%            are cutpnt(i-1) < y <= cutpnt(i)
%    .hot0 : list of residuals greater than 3.25 standard deviations from the median
%      .ho : list of the values of residuals greater than 3.25
%      .me : mean of e
%    .rstd : standard deviation of the mean of e
%   .rtval : t-value of the mean of e
%    .maxe : maximum value of e
%     .mde : median value of e
%    .mine : minimum value of e
%    .skew : skewness
%    .kurt : kurtosis
%     .bst : Bowman-Shenton normality statistic
%     .pnt : p-value of the Bowman-Shenton statistic
%     .tsk : p-value of skewness
%     .tkr : p-value of kurtosis
%      .dw : Durbin-Watson statistic
%     .tdw : t-value of the Durbin-Watson statistic
%    .ptdw : p-value of the Durbin-Watson statistic
%      .n0 : number of residuals lower than the median
%      .n1 : number of residuals higher than the median
%      .nr : number of runs on residuals
%    .Tval : t-value of the number of runs on residuals
%      .rs : autocorrelations of squared residuals
%     .pcs : partial correlations of squared residuals
%  .qstats : Q-statistics based on squared residuals
%   .pvals : p-values of the Q-statistics based on squared residuals
%     .dfs : degrees of freedom for the Q-statistics based on squared
%            residuals
%    .seas : standard errors associated with squared residuals
%       .h : closest integer to ne/3 needed in the computation of H;
%            degrees of freedom of the F-distribution
%       .H : heteroskedasticity statistic
%      .pH : p-value of the heteroskedasticity statistic
%     .aic : Akaike information criterion
%     .bic : Bayes information criterion
%    .mser : mean squared error of residuals
%   .stder : standard error of residuals
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


dn = double(ndrs); %information criteria
dnp = double(nr+nreg);
aic = dn * (log(2*pi) + log(Ff)) + 2 * dnp;
bic = dn * (log(2*pi) + log(Ff)) + log(dn) * dnp;

ne = length(e); %residual sample length
ve = var(e, 1);
stde = sqrt(ve); %residual variance and std. dev.
%the following line added 21-11-2009
mser = Ss / (dn - nr - nreg);
stder = sqrt(ve); %mse and standard error of resid
me = mean(e); %residual mean
%the follwing line commented 21-11-2009
% ec=e-me*ones(size(e));                           %residuals must be mean-corrected before autocv
%the following line added 21-11-2009
ec = e;
ic = 1;
[c0, cv, r] = autcov(ec, lag, ic); %autocorrelations
[fi, pc] = durlev(c0, cv); %partial autocorrelations
orders = 1:lag;
[qstat, pval, df, sea] = lbs(ne, orders, r, nr);
sep = ones(lag, 1) / sqrt(ne);
rstd = sqrt(ve/ne);
rtval = me / rstd; %standard deviation and t-value of mean
maxe = max(e);
mde = median(e);
mine = min(e); %median, minimum and maximum
[skew, kurt, sk, ss] = skewkur(e, me, ve, 0, ne); %skweness and kurtosis
dw = durwat(e, 0, ne, Ss); %Durbin-Watson
tsk = skew / ss;
tkr = (kurt - 3) / sk;
bst = tsk^2 + tkr^2; %normality test
sqe = e.^2;
h = round(ne/3);
H = sum(sqe(ne-h+1:ne)) / ...
    sum(sqe(1:h));%heteroscedasticity test (Harvey, 1987, p. 259)
pnt = 1 - gammp(double(2)*.5d0, bst*.5d0); %p-values
tsk = 1 - gammp(double(1)*.5d0, (tsk^2)*.5d0);
tkr = 1 - gammp(double(1)*.5d0, (tkr^2)*.5d0);
pf = fdis_cdf(H, h, h);
if (pf < .5)
    pH = 2 * pf;
else
    pH = 2 * (1 - pf);
end
ddw = 4.d0 / double(ne);
tdw = (dw - 2.d0) / sqrt(ddw);
[Result, Ccum] = cumnor(tdw);
if (tdw < 0)
    ptdw = 2 * Result;
else
    ptdw = 2 * Ccum;
end
[n0, n1, nr, Tval] = runcom(e, ne, mde);

[bin, cutpnt, otlrt0, otlr] = hist2(e, mde); %residual histogram
no = bin;
xo = cutpnt;

sqe = sqe - mean(sqe) * ones(size(sqe)); %squared residuals
[c0s, cvs, rs] = autcov(sqe, lag, ic); %autocorrelations of squared residuals
[fi, pcs] = durlev(c0s, cvs); %partial autocorrelations of squared residuals
[qstats, pvals, dfs, seas] = lbs(ne, orders, rs, 0);


%end of residual computations


%put all residual results in a structure
infr.e = e;
infr.ne = ne;
infr.ve = ve;
infr.stde = stde;
infr.sconp = sconp;
infr.conp = conp;
infr.orders = orders;
infr.r = r;
infr.pc = pc;
infr.qstat = qstat;
infr.pval = pval;
infr.df = df;
infr.sea = sea;
infr.sep = sep;
infr.no = no;
infr.xo = xo;
infr.hot0 = otlrt0;
infr.ho = otlr;
infr.me = me;
infr.rstd = rstd;
infr.rtval = rtval;
infr.maxe = maxe;
infr.mde = mde;
infr.mine = mine;
infr.skew = skew;
infr.kurt = kurt;
infr.dw = dw;
infr.bst = bst;
infr.pnt = pnt;
infr.tsk = tsk;
infr.tkr = tkr;
infr.dw = dw;
infr.tdw = tdw;
infr.ptdw = ptdw;
infr.n0 = n0;
infr.n1 = n1;
infr.nr = nr;
infr.Tval = Tval;
infr.rs = rs;
infr.pcs = pcs;
infr.qstats = qstats;
infr.pvals = pvals;
infr.dfs = dfs;
infr.seas = seas;
infr.h = h;
infr.H = H;
infr.pH = pH;
infr.aic = aic;
infr.bic = bic;
infr.mser = mser;
infr.stder = stder;
%end of residual structure
