function mprintr(result, fid)
% PURPOSE: print a two column table with the estimates and their t-values
%          contained in structure result
%---------------------------------------------------
% USAGE:     mprint(x,info)
% where: fid       = file identifier for output (default = 1)
%        result    = a structure containing estimation results with the
%                    following fields:
%           .xvf : estimated parameters
%           .xf : vector of fixed parameters
%      .sigma2c : concentrated parameter estimate
%       .Sigmar : estimated exact covariance matrix of residuals
%           .tv : t-values of the estimated varma parameters
%    .residexct : matrix containing recursive residuals, only if Y is empty
%            .e : vector of standardized residuals at the end of estimation
%                 (Q'_2*y)
%           .ff : vector of nonlinear functions whose sum of squares is
%                 minimized at the end of estimation
%           .h  : vector of estimated regression estimates
%           .H  : matrix of mse of h
%           .A  : estimated state vector, x_{t|t-1}, obtained with the
%                 Kalman filter at the end of the sample
%           .P  : Mse of A
%           .tvr: vector of t-values for h
%       .ferror : flag for errors
%---------------------------------------------------
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
%

if (nargin < 2)
    fid = 1;
end
z = [result.xvf', result.tv];
clear in
in.cnames = char('  Estimate', '   T-ratio');
m = length(result.xvf);
tit = 'Parameter           ';
for i = 1:m;
    tit = char(tit, ['par. #', int2str(i)]);
end
in.rnames = tit;
in.fmt = char('%12.4f', '%12.4f');
in.fid = fid;
mprint(z, in);
if ~isempty(result.h)
    z = [result.h, result.tvr];
    clear in
    fid = 1;
    in.cnames = char('  Estimate', '   T-ratio');
    m = length(result.h);
    tit = 'Regression Parameter';
    for i = 1:m;
        tit = char(tit, ['par. #', int2str(i)]);
    end
    in.rnames = tit;
    in.fmt = char('%12.4f', '%12.4f');
    in.fid = fid;
    disp(' ');
    mprint(z, in);
end
