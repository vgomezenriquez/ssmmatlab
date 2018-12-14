function outa = usmestos(fname, fmeta)
%
% Function for automatic identification, estimation and forecasting of
% a structural model for one series
%
%     INPUTS:
%      fname : If fmeta = 0 or absent, a string such that fname.m is a
%              matlab function in the spec subdirectory that returns the
%              structure ser. In this structure, instruction for this
%              function are given. If fmeta = 1, a string such that
%              fname.txt contains a list of names of matlab functions in
%              the spec subdirectory that will be treated sequentially.
%      fmeta : = 0, fname.m is a matlab function in the spec subdirectory
%              that returns the structure ser; = 1, fname.txt is a file
%              that contains a list of matlab functions in the spec
%              subdirectory that will be treated sequentially. If not
%              input, the program sets by default fmeta = 0,
%
%  OUTPUTS:
%      outa  : a structure containing model information for the input
%              with fields:
%       title: a string with the name of series
%      nziyip: a 1 x 3 array with number of obs., initial year, initial per.
%        freq: number of seasons
%        orig: original series
%       model: structre with model information. It contains the following
%              fields
%           lam: flag for logarithmic transformation, = 0, take logs, = 1,
%                do not take logs
%             X: X matrix in the state space form
%             Z: Z matrix in the state space form
%             G: G matrix in the state space form
%             W: W matrix in the state space form
%             T: T matrix in the state space form
%             H: H matrix in the state space form
%           ins: ins matrix for the initial conditions
%             i: i array for the initial conditions
%        resinf: structure containing residual information
%         sconp: residual standard error
%       StochCc: matrix containing the stochastic components
%      StochSCc: matrix containing the mse of the stochastic components
%      oStochCc: matrix containing the stochastic components in the
%                original scale
%     oStochSCc: matrix containing the mse of the stochastic components in
%                the original scale
%            Cc: matrix containing the components including deterministic
%                effects
%           SCc: matrix containing the mse of Cc
%           oCc: matrix containing the Cc in the original scale
%          oSCc: matrix containing the mse of the oCc
%           npr: number of forecasts
%            Xp: matrix containing the forecasts of X
%            Wp: matrix containing the forecasts of W
%           pry: forecasts
%          spry: mse of the forecasts
%          alpr: matrix containing the forecasts of the state vector
%         malpr: three dimensional array containing each of the covariance
%                matrices of alpr
%         salpr: matrix containing the mse of alpr
%          opry: forecasts in the original scale
%         ospry: mse of the forecasts in the original scale
%         oalpr: matrix containing the alpr in the original scale
%        osalpr: matrix containing the mse of oalpr
%         ser: the input structure
%      result: a structure containing estimation results. It has
%              the following fields:
%          xvf: array containing the estimated parameters
%           xf: array containing the fixed parameters
%            e: array containing the residuals
%           Ss: residual sum of squares
%           Ff: the product F'*F
%      sigma2c: standard error of the parameter concentrated out of the
%               likelihood
%         Pevf: prediction error variance
%        SPevf: square root of Pevf
%           tv: t-values of the estimated parameters
%           se: standard errors of the estimated parameters
%           .F: vector of nonlinear functions whose sum of squares is
%               minimized at the end of estimation
%            h: vector of regression parameters
%            M: mse of h
%            A: Augmented state vector at the end of filtering
%            P: mse matrix of A at the end of filtering
%       olsres: OLS residuals
%          tvr: t-values of the regression parameters
%          ser: standard errors of the regression parameters
%       ferror: flag for errors
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

% tic
% profile on       %this is to optimize the code. It gives some information
%about the performance of all of the functions.

if nargin < 2
    fmeta = 0; %flag for an external file
end

if fmeta == 1
    mname = ['usmspec', filesep, fname, '.txt'];
    fid = fopen(mname, 'r');
    nmatr = fgetl(fid);
    while feof(fid) == 0
        tline = fgetl(fid);
        nmatr = char(nmatr, tline);
    end
else
    nmatr = fname;
end
% nmatr

outa = [];

N = size(nmatr, 1);
%file for summary report
% frname=fullfile('results','summary.txt');
% fidr = fopen(frname,'w');

pathc = pwd; %current directory
addpath([pathc, filesep, 'usmspec']) %add specusm subdirectory to the path

for ii = 1:N
    fname = nmatr(ii, :);
    fprintf(1, '%s\n', fname)
    
    dbname = deblank(fname);
    clear ser %structure containing series information
    
    %read specifications and data from external files
    
    ser = feval(dbname);
    
    outa = [outa, usmestni(dbname, ser)];
    
end

% fclose(fidr); %close summary file
% profile viewer
% profile off
% toc
