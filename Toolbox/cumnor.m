function [Result, Ccum] = cumnor(Arg)
%     Last change:  BCM  21 Nov 97   10:07 pm
%      SUBROUTINE cumnor(Arg,Result,Ccum)
%      IMPLICIT NONE
%C**********************************************************************
%C
%C     SUBROUINE CUMNOR(X,RESULT,CCUM)
%C
%C
%C                              Function
%C
%C
%C     Computes the cumulative  of    the  normal   distribution,   i.e.,
%C     the integral from -infinity to x of
%C          (1/sqrt(2*pi)) exp(-u*u/2) du
%C
%C     X --> Upper limit of integration.
%C                                        X is DOUBLE PRECISION
%C
%C     RESULT <-- Cumulative normal distribution.
%C                                        RESULT is DOUBLE PRECISION
%C
%C     CCUM <-- Compliment of Cumulative normal distribution.
%C                                        CCUM is DOUBLE PRECISION
%C
%C
%C     Renaming of function ANORM from:
%C
%C     Cody, W.D. (1993). "ALGORITHM 715: SPECFUN - A Portabel FORTRAN
%C     Package of Special Function Routines and Test Drivers"
%C     acm Transactions on Mathematical Software. 19, 22-32.
%C
%C     with slight modifications to return ccum and to deal with
%C     machine constants.
%C
%C**********************************************************************
%C
%C
%C Original Comments:
%C------------------------------------------------------------------
%C
%C This function evaluates the normal distribution function:
%C
%C                              / x
%C                     1       |       -t*t/2
%C          P(x) = ----------- |      e       dt
%C                 sqrt(2 pi)  |
%C                             /-oo
%C
%C   The main computation evaluates near-minimax approximations
%C   derived from those in "Rational Chebyshev approximations for
%C   the error function" by W. J. Cody, Math. Comp., 1969, 631-637.
%C   This transportable program uses rational functions that
%C   theoretically approximate the normal distribution function to
%C   at least 18 significant decimal digits.  The accuracy achieved
%C   depends on the arithmetic system, the compiler, the intrinsic
%C   functions, and proper selection of the machine-dependent
%C   constants.
%C
%C*******************************************************************
%C*******************************************************************
%C
%C Explanation of machine-dependent constants.
%C
%C   MIN   = smallest machine representable number.
%C
%C   EPS   = argument below which anorm(x) may be represented by
%C           0.5  and above which  x*x  will not underflow.
%C           A conservative value is the largest machine number X
%C           such that   1.0 + X = 1.0   to machine precision.
%C*******************************************************************
%C*******************************************************************
%C
%C Error returns
%C
%C  The program returns  ANORM = 0     for  ARG .LE. XLOW.
%C
%C
%C Intrinsic functions required are:
%C
%C     ABS, AINT, EXP
%C
%C
%C  Author: W. J. Cody
%C          Mathematics and Computer Science Division
%C          Argonne National Laboratory
%C          Argonne, IL 60439
%C
%C  Latest modification: March 15, 1992
%C
%C------------------------------------------------------------------
%      INTEGER i
%      DOUBLE PRECISION a,Arg,b,c,d,del,eps,half,p,one,q,Result,sixten,
%     &                 temp,sqrpi,thrsh,root32,x,xden,xnum,y,xsq,zero,
%     &                 minx,Ccum
%      DIMENSION a(5),b(4),c(9),d(8),p(6),q(5)
%C------------------------------------------------------------------
%C  External Function
%C------------------------------------------------------------------
%      DOUBLE PRECISION spmpar
%      EXTERNAL spmpar
%C------------------------------------------------------------------
%C  Mathematical constants
%C
%C  SQRPI = 1 / sqrt(2*pi), ROOT32 = sqrt(32), and
%C  THRSH is the argument for which anorm = 0.75.
%C------------------------------------------------------------------
one = 1.0d0;
half = 0.5d0;
zero = 0.0d0;
sixten = 1.60d0;
sqrpi = 3.9894228040143267794d-1;
thrsh = 0.66291d0;
root32 = 5.656854248d0;
%C------------------------------------------------------------------
%C  Coefficients for approximation in first interval
%C------------------------------------------------------------------
a = [2.2352520354606839287d00, 1.6102823106855587881d02, ...
    1.0676894854603709582d03, 1.8154981253343561249d04, ...
    6.5682337918207449113d-2];
b = [4.7202581904688241870d01, 9.7609855173777669322d02, ...
    1.0260932208618978205d04, 4.5507789335026729956d04];
%C------------------------------------------------------------------
%C  Coefficients for approximation in second interval
%C------------------------------------------------------------------
c = [3.9894151208813466764d-1, 8.8831497943883759412d00, ...
    9.3506656132177855979d01, 5.9727027639480026226d02, ...
    2.4945375852903726711d03, 6.8481904505362823326d03, ...
    1.1602651437647350124d04, 9.8427148383839780218d03, ...
    1.0765576773720192317d-8];
d = [2.2266688044328115691d01, 2.3538790178262499861d02, ...
    1.5193775994075548050d03, 6.4855582982667607550d03, ...
    1.8615571640885098091d04, 3.4900952721145977266d04, ...
    3.8912003286093271411d04, 1.9685429676859990727d04];
%C------------------------------------------------------------------
%C  Coefficients for approximation in third interval
%C------------------------------------------------------------------
p = [2.1589853405795699d-1, 1.274011611602473639d-1, ...
    2.2235277870649807d-2, 1.421619193227893466d-3, ...
    2.9112874951168792d-5, 2.307344176494017303d-2];
q = [1.28426009614491121d00, 4.68238212480865118d-1, ...
    6.59881378689285515d-2, 3.78239633202758244d-3, ...
    7.29751555083966205d-5];
%C------------------------------------------------------------------
%C  Machine dependent constants
%C------------------------------------------------------------------
%     eps=spmpar(1)*0.5d0;
minx = realmin; %minx=spmpar(2);
%C------------------------------------------------------------------
x = Arg;
y = abs(x);
if (y < thrsh)
    %C------------------------------------------------------------------
    %C  Evaluate  anorm  for  |X| <= 0.66291
    %C------------------------------------------------------------------
    xsq = zero;
    if (y > eps), xsq = x * x;
    end
    xnum = a(5) * xsq;
    xden = xsq;
    for i = 1:3
        xnum = (xnum + a(i)) * xsq;
        xden = (xden + b(i)) * xsq;
    end
    Result = x * (xnum + a(4)) / (xden + b(4));
    temp = Result;
    Result = half + temp;
    Ccum = half - temp;
    %C------------------------------------------------------------------
    %C  Evaluate  anorm  for 0.66291 <= |X| <= sqrt(32)
    %C------------------------------------------------------------------
elseif (y <= root32)
    xnum = c(9) * y;
    xden = y;
    for i = 1:7
        xnum = (xnum + c(i)) * y;
        xden = (xden + d(i)) * y;
    end
    Result = (xnum + c(8)) / (xden + d(8));
    xsq = fix(y*sixten) / sixten; %xsq=aint(y*sixten)/sixten;
    del = (y - xsq) * (y + xsq);
    Result = exp(-xsq*xsq*half) * exp(-del*half) * Result;
    Ccum = one - Result;
    if (x > zero)
        temp = Result;
        Result = Ccum;
        Ccum = temp;
    end
    %C------------------------------------------------------------------
    %C  Evaluate  anorm  for |X| > sqrt(32)
    %C------------------------------------------------------------------
else
    Result = zero;
    xsq = one / (x * x);
    xnum = p(6) * xsq;
    xden = xsq;
    for i = 1:4
        xnum = (xnum + p(i)) * xsq;
        xden = (xden + q(i)) * xsq;
    end
    Result = xsq * (xnum + p(5)) / (xden + q(5));
    Result = (sqrpi - Result) / y;
    xsq = fix(x*sixten) / sixten; %xsq=aint(x*sixten)/sixten;
    del = (x - xsq) * (x + xsq);
    Result = exp(-xsq*xsq*half) * exp(-del*half) * Result;
    Ccum = one - Result;
    if (x > zero)
        temp = Result;
        Result = Ccum;
        Ccum = temp;
    end
end
if (Result < minx), Result = 0.0d0;
end
if (Ccum < minx), Ccum = 0.0d0;
end
%C------------------------------------------------------------------
%C  Fix up for negative argument, erf, etc.
%C------------------------------------------------------------------
%C----------Last card of ANORM ----------
%     end
