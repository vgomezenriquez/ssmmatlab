function infm = minfm(f, tr, mvx, tolf, maxit, nu0, jac, prt, chb, inc)
%*************************************************************************
% This function puts the function name and optimization options into
% structure infm
%
%   INPUTS:
%       f : a function to evaluate the vector ff of individual functions
%           such that ff'*ff is minimized
%      tr > 0 : x is passed from marqdt to f but not passed from f to marqdt
%         = 0 : x is passed from marqdt to f and passed from f to marqdt
%     mvx :  =1 exact maximum likelihood
%            =0 unconditional least squares
%    tolf : parameter used for stopping
%   maxit : maximum number of iterations
%     nu0 : initial value of the nu parameter
%     jac = 1 : evaluation of jacobian and gradient at the solution is performed
%         = 0 : no evaluation of jacobian and gradient at the solution is performed
%     prt = 1 : printing of results
%         = 0 : no printing of results
%     chb = 1 : compute the beta estimate and its MSE
%           0 : do not compute the beta estimate and its MSE
%     inc = 0, the initial states in the filter equations to obtain the
%              filtered variables are equal to zero (not estimated)
%         = 1, the initial states in the filter equations are estimated
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


infm.f = f;
infm.tr = tr;
infm.mvx = mvx;
infm.tolf = tolf;
infm.maxit = maxit;
infm.nu0 = nu0;
infm.jac = jac;
infm.prt = prt;
infm.chb = chb;
infm.inc = inc;
