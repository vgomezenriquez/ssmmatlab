function [Skew, Kurt, Sk, Ss] = skewkur(e, me, ve, nr, ne)
%**************************************************************************
% This function calculates components of the skewness and kurtosis test
% statistics
%
%     INPUTS:
%         e : residual vector
%        me : mean of e
%        ve : variance of e
%        nr : number of regression variables
%        ne : length of e
%
%    OUTPUTS:
%      Skew : skewness
%      Kurt : kurtosis
%        Sk : constant to obtain kurtosis test statistic: (Kurt/Sk)^2
%        Ss : constant to obtain skewness test statistic: (Skew/Ss)^2
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
%**************************************************************************


den = double(ne-nr);
% Skew=0; Kurt=0;
% for i=1+nr:ne
%  Skew=Skew+(e(i)-me)^3;
%  Kurt=Kurt+(e(i)-me)^4;
% end
Skew = sum((e(1+nr:ne) - me * ones(size(e(1+nr:ne)))).^3);
Kurt = sum((e(1+nr:ne) - me * ones(size(e(1+nr:ne)))).^4);
Skew = Skew / (ve^1.5 * den);
Kurt = Kurt / (ve^2 * den);
Ss = sqrt(6.d0/den);
Sk = sqrt(24.d0/den);
