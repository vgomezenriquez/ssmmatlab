function [x, J] = arimaopt(fmarqdt, fid, x0, xv, xf, y, Y, parm, infm, pr)
%
% This function performs the optimization for an ARIMA model
%
%     INPUTS:
%     fmarqdt : = 1 estimation with lsqnonlin (matlab), = 0, estimation
%               with marqdt
%     fid     : an integer, corresponding to an external file
%     x0      : an array containing the initial estimated values
%     xv      : an array containing the variable parameters
%     xf      : an array containing the fixed parameters
%     y       : input series
%     Y       : matrix of regression variables
%     parm: a structure containing model information, where
%       .s:  seasonality
%       .S:  second seasonality
%       .p:  AR order
%      .ps: order of the AR of order s
%       .q:  order of the regular MA
%      .qs: order of the MA of order s (1 at most)
%      .qS: order of the MA of order S (1 at most)
%      .dr: order of regular differencing
%      .ds: order of differencing of order s
%      .dS: order of differencing of order S
%    .pvar:  array containing the indices of variable parameters
%    .pfix:  array containing the indices of fixed parameters
%  .ninput: number of inputs
%  .inputv: array containing the lagged input variables corresponding to
%           the polynomial approximations to the rational input filters
%   .delay: array with the delays of the input filters
%      .ma: array with the ma parameters of the input filters
%      .ar: array with the ar parameters of the input filters
%     .npr: number of forecasts
%   infm      : structure containing function names and optimization options
%   .f  :   a function to evaluate the vector ff of individual functions
%           such that ff'*ff is minimized
%   .tr :   >0 x is passed from marqdt to f but not passed from f to marqdt
%           =0 x is passed from marqdt to f and passed from f to marqdt
%   .tol:   a parameter used for stopping
%   .jac:   =1 evaluation of jacobian and gradient at the solution is performed
%           =0 no evaluation of jacobian and gradient at the solution is performed
% .maxit:   maximum number of iterations
%   .nu0:   initial value of the nu parameter
%   .prt:   =1 printing of results
%           =0 no printing of results
%   .chb:   = 1  compute the beta estimate and its MSE
%             0  do not compute the beta estimate and its MSE
%   .inc:   = 0, the initial states in the filter equations to obtain the
%                filtered variables are equal to zero (not estimated)
%           = 1, the initial states in the filter equations are estimated
%     pr      : = 1, print results in an external file, = 0, do not print
%
%  OUTPUTS:
%      x     : an array containing the estimated parameter values
%      J     : a matrix containing the Jacobian at the solution
%
%
% Copyright (c) 21 July 2015 by Victor Gomez
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

if fmarqdt == 0
    if pr == 1, fprintf(fid, 'Arima Estimation (Lsqnonlin):\n');
    end
    %     tic;
    [x, J, exitflag] = optim1(xv, y, Y, parm, infm);
    %     time = toc;
    %     fprintf(fid,'Elapsed time: %8.2f\n\n',time);
else
    %     tic;
    [x, J, ff, g, iter, conf] = marqdt(infm, xv, y, Y, parm, infm, xf);
    %     time = toc;
    if pr == 1
        fprintf(fid, 'Arima Estimation (Levenberg-Marquardt):\n'); %Levenberg-Marquardt
        fprintf(fid, 'Number of iterations: %5i\n', iter);
        fprintf(fid, 'Number of function evaluations: %5i\n', conf);
        %     fprintf(fid,'Elapsed time: %8.2f\n\n',time);
    end
end
