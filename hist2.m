function [bin, cutpnt, otlrt0, otlr] = hist2(Y, med)
%
%
% c-----------------------------------------------------------------------
% c     Calculates the histogram of the Nobs long data vector, y.
% c-----------------------------------------------------------------------
%
% after calling this function:
%  bar(cutpnt,bin)
% xlabel('standard deviation intervals')
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
% c bin     i  Local vector of counts of observations between cut points
% c             i.e. bins
% c cutpnt  d  Local vector of cut points where observations counted in
% c             bin(i) are cutpnt(i-1) < y <= cutpnt(i)
% c             nobs-nefobs
% c i       i  Local do loop index
% c ibin    i  Local index for the current bin, or column of the
% c             histogram
% c irow    i  Local index for the current row of the output
% c lowbnd  d  Local scalar for the low bound of the observations
% c median    d  Local median of the y's
% c nbin    i  Local number of bins or columns in the histogram
% c notlr   i  Local number of outliers (in otlr)
% c nobs    i  Input number of observations
% c nrow    i  Local number of rows in the histogram
% c otlr    i  Local notlr long list of the values of residuals
% c             greater than 3.25
% c otlrt0  i  Local notlr long list of residuals greater than 3.25
% c             standard deviations from the median
% c scale   d  Local scale factor to make the number of rows be 40
% c stddev  d  Local standard deviation of the y's
% c tmp     d  Local temporary scalar
% c width   d  Local width between cut points
% c y       d  Input nobs long vector of observations
% c-----------------------------------------------------------------------
% c     ------------------------------------------------------------------
% c-----------------------------------------------------------------------
% c     Find the minimum, maximum, median, and standard deviation of the
% c observations.
% c-----------------------------------------------------------------------
%   Mt1=7;

if nargin < 2
    med = median(Y);
end
Nobs = length(Y);
%  srtdy=abs(Y);
%  srtdy=median(srtdy);
% C
% C MAD ESTIMATOR OF RESIDUAL STANDARD DEVIATION
% C
stddev = 1.49D0 * std(Y, 1);
% c-----------------------------------------------------------------------
% c     Find the range and lower bound of the histogram.
% c-----------------------------------------------------------------------
% C     width=.5D0
% C     lowbnd=-3.25D0
% C     nbin=15
if Nobs <= 100
    lowbnd = -3.D0;
    nbin = 13;
elseif Nobs <= 250
    lowbnd = -3.5D0;
    nbin = 15;
else
    lowbnd = -4.D0;
    nbin = 17;
end
width = 2.D0 * (-lowbnd) / double(nbin-2);
cutpnt = zeros(nbin, 1);
otlr = zeros(ceil(Nobs/4));
otlrt0 = zeros(ceil(Nobs/4));
% c-----------------------------------------------------------------------
% c     Calculate the cut points.
% c-----------------------------------------------------------------------
tmp = lowbnd;
% c     ------------------------------------------------------------------
for i = 1:nbin
    cutpnt(i) = tmp;
    tmp = tmp + width;
end
% c-----------------------------------------------------------------------
% c     Sort the observations into bins.
% c-----------------------------------------------------------------------
nrow = 0;
notlr = 0;
bin = zeros(nbin, 1);
%  c     ------------------------------------------------------------------
for i = 1:Nobs
    flag = 0;
    tmp = (Y(i) - med) / stddev;
    % c     ------------------------------------------------------------------
    if tmp < lowbnd
        notlr = notlr + 1;
        otlrt0(notlr) = i;
        otlr(notlr) = tmp;
        bin(1) = bin(1) + 1;
        nrow = max(nrow, bin(1));
        % c     ------------------------------------------------------------------
    else
        for ibin = 2:nbin - 1
            if tmp < cutpnt(ibin)
                bin(ibin) = bin(ibin) + 1;
                nrow = max(nrow, bin(ibin));
                flag = 1;
                break
            end
        end
        if flag == 0
            % c     ------------------------------------------------------------------
            notlr = notlr + 1;
            otlrt0(notlr) = i;
            otlr(notlr) = tmp;
            bin(nbin) = bin(nbin) + 1;
            nrow = max(nrow, bin(nbin));
        end
    end
end
% c     ------------------------------------------------------------------
scale = ceil(double(nrow)/100.D0);
if scale > 1
    for ibin = 1:nbin
        bin(ibin) = floor(bin(ibin)/scale);
    end
end
if notlr == 0
    otlr = [];
    otlrt0 = [];
else
    otlr = otlr(1:notlr);
    otlrt0 = otlrt0(1:notlr);
end
% c-----------------------------------------------------------------------
% c     Print the histogram sideways with negative values on top.
% c-----------------------------------------------------------------------
%       WRITE(Mt1,1010)
%  1010 FORMAT(//'  RESIDUAL HISTOGRAM ',/
%      &'  STANDARD',/,'  DEVIATION',/,
%      &'  INTERVALS      FREQUENCY')
% c     ------------------------------------------------------------------
%       DO irow=1,nbin-1
%         WRITE(Mt1,1030)cutpnt(irow),('X',i=1,bin(irow))
%  1030   FORMAT('   < ',F5.2,' |',69A1)
%       END DO
% c     ------------------------------------------------------------------
%   WRITE(Mt1,1090)cutpnt(nbin-1),('X',i=1,bin(nbin))
%  1090 FORMAT('   > ',F5.2,' |',69A1)
% c     ------------------------------------------------------------------
%       WRITE(Mt1,1050)int(scale)
%  1050 FORMAT(/,'  ONE ''X''=',I2,' OBSERVATION(S)')
% c     ------------------------------------------------------------------
%       IF(notlr.gt.0)THEN
%        WRITE(Mt1,1060)cutpnt(nbin-1)
%  1060  FORMAT(/,'  RESIDUALS WITH |T|>',F5.2,/,'  OBS       T-VALUE',/,
%      &        '  -----------------')
%        DO i=1,notlr
%         WRITE(Mt1,1070)otlrt0(i),otlr(i)
%  1070   FORMAT(' ',I4,T12,F8.2)
%        END DO
%       END IF
% c     ------------------------------------------------------------------
%       WRITE(Mt1,1080)med,stddev
%  1080 FORMAT(/,'  MEDIAN',T25,F15.3,/,'  ROBUST STD. DEV.(MAD)',
%      &       T25,F15.3)
% c     ------------------------------------------------------------------
