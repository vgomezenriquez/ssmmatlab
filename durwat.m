function Dw = durwat(Res, J, K, Ss)
%**************************************************************************
% This function computes the Durbin-Watson statistic
%
%   INPUTS:
%     Res : residual vector
%       J : number of first residuals in Res not used in the
%           computation of Dw
%       K : integer specifying the last residual in Res used in the
%           computation of Dw
%      Ss : variance of Res
%
%   OUTPUT:
%      Dw : Durbin-Watson statistic
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

Dw = sum((Res(J+2:K) - Res(J+1:K-1)).^2) / Ss;
