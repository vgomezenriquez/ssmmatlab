function outa = arimaestos(fname, fmeta)
%
% Function for automatic identification, estimation and forecasting of
% ARIMA or transfer function models for one or several series
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
%       model: structre with ARIMA model information. In the case of a
%              transfer function, it is the ARIMA model corresponding to
%              the finite linear approximation to the input filters. It has
%              the following fields:
%          lam: = 0, logs are taken, = 1, no logs
%         mean: = 1, a mean is added, = 0, no mean
%            p: degree of regular AR polynomial
%            d: degree of regular differencing
%            q: degree of regular MA polynomial
%           ps: degree of seasonal AR polynomial
%           ds: degree of seasonal differencing
%           qs: degree of seasonal MA polynomial
%         nreg: the number of regression variables
%       result: a structure containing estimation results
%          phi: an array containing the regular AR polynomial coefficients
%         phis: an array containing the seasonal AR polynomial coefficients
%           th: an array containing the regular MA polynomial coefficients
%          ths: an array containing the seasonal MA polynomial coefficients
%        nrout: number of outliers
%            C: critical value for outlier detection
%         nind: observation numbers of the outliers
%          tip: string containing the outlier types
%       matsis: a structure containing the state space form of the model
%       resinf: a structure containing information about the residuals
%           hb: array containing the regression estimates
%           Mb: matrix containing the covariance matrix of the regression
%               estimates
%            Y: matrix containing the total regression effects
%          seb: array containing the regression standard errors
%           tb: array containing the t-values of the regression estimates
%          Yrg: array containing the regression variables that are not
%               outliers
%        Youtg: array containing the outlier variables
%           se: array containing the standard errors of the estimates
%           tt: array containing the t-values of the estimates
%          npr: number of forecasts
%          pry: array containing the forecasts (transformed scale)
%         spry: array containing the standard errors of the forecasts
%         opry: same as pry but in the original scale
%        ospry: same os spry but in the original scale
%     tfmodel: structure with transfer function model information.
%       matsis: a structure containing the state space form of the model
%               correponding to the filtered inputs
%       result: a structure containing estimation results
%         nreg: the number of regression variables
%          Yrg: array containing the regression variables that are not
%               outliers
%        Youtg: array containing the outlier variables
%          yci: output corrected by filtered inputs
%        tford: a three column array in which the i-th row has three
%               numbers corresponding to the delay, the numerator degree
%               and the denominator degree of the i-th input filter
%          phi: an array containing the regular AR polynomial coefficients
%         phis: an array containing the seasonal AR polynomial coefficients
%           th: an array containing the regular MA polynomial coefficients
%          ths: an array containing the seasonal MA polynomial coefficients
%          omg: a cell array containing the numerators of the input filters
%          del: a cell array containing the denominators of the input
%               filters
%       resinf: a structure containing information about the residuals
%           se: array containing the standard errors of the estimates
%           tt: array containing the t-values of the estimates
%          npr: number of forecasts
%         dpry: array containing the forecasts of the output corrected by
%               the filtered inputs
%        dspry: array containing the standard errors of the forecasts of
%               the output corrected by the filtered inputs
%          Yin: array containing the input variables
%      modpred: a multiple structure containing the input forecasts if any.
%               the forecasts for each input are given in field .pred
%     modinput: a multiple structure containing the models for the inputs
%               if any (.mod = 0, no model; .mod =1, there is model). The
%               model for each input has fields .alpha, .phi, .theta,
%               .sigma2
%            y: output series in the transformed scale
%          pry: array containing the forecasts (transformed scale)
%         spry: array containing the standard errors of the forecasts
%         opry: same as pry but in the original scale
%        ospry: same os spry but in the original scale
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

% tic
% profile on       %this is to optimize the code. It gives some information
%about the performance of all of the functions.

if nargin < 2
    fmeta = 0; %flag for an external file
end

if fmeta == 1
    mname = ['spec', filesep, fname, '.txt'];
    fid = fopen(mname, 'r');
    nmatr = [];
    while feof(fid) == 0
        tline = fgetl(fid);
        nmatr = strvcat(nmatr, tline);
    end
else
    nmatr = fname;
end
% nmatr

outa = [];

[N, M] = size(nmatr);
%file for summary report
frname = fullfile('results', 'summary.txt');
fidr = fopen(frname, 'w');

pathc = pwd; %current directory
addpath([pathc, filesep, 'spec']) %add spec subdirectory to the path

for ii = 1:N
    fname = nmatr(ii, :);
    fprintf(1, '%s\n', fname)
    
    dbname = deblank(fname);
    clear ser %structure containing series information
    
    %read specifications and data from external files
    
    ser = feval(dbname);
    
    if ~isfield(ser, 'ninput')
        outa = [outa, arimaestni(dbname, ser, fidr, ii)];
    else
        outa = [outa, arimaestwi(dbname, ser, fidr, ii)];
    end
    
end

fclose(fidr); %close summary file
% profile viewer
% profile off
% toc
