function [fjac, g] = fdjac2(info, x, ff, varargin)
%
% This function computes a finite-difference jacobian
%
% Input arguments:
%   info structure containing function names and optimization options
%   .f  :   a function to evaluate the vector ff of individual functions
%           such that ff'*ff is minimized
%   .tr :   >0 a transformation of the parameters corresponding to tr variables
%              is performed
%           =0 no transformation of parameters is performed
%   .tol:   a parameter used for stopping
%   .jac:   =1 evaluation of jacobian at gradient at the solution is performed
%           =0 no evaluation of jacobian at gradient at the solution is performed
%   .max:   maximum number of iterations
%   .nu0:   initial value of the nu parameter
%   .prt:   =1 printing of results
%           =0 no printing of results
% ff:       the vector of functions evaluated at x
% x:        a vector containing the initial parameter values
% varargin: arguments to be passed to function f
%
% Output arguments:
% fjac     a matrix containing the jacobian
% g        a vector containing the (1/2)gradient
%

eps = 1.19 * 1e-7; %significant digits

% h = eps^(1/3)*max(abs(x),1);
h = sqrt(eps) * max(.1, abs(x)); %stepsize modified by Victor Gomez, June, 2003
for i = 1:length(h)
    if x(i) < 0
        h(i) = -h(i);
    end
end
n = length(x);
k = length(ff);
fjac = zeros(k, n);
fHandle = str2func(info.f);
%  for cellfun compute
for j = 1:n;
    x1 = x;
    x1(j) = x(j) + h(j);
    %  1 dec 2014
    
    [ffd] = fHandle(x1, varargin{:});
    %      ffd=feval(info.f,x1,varargin{:});
    fjac(:, j) = (ffd - ff) / h(j);
end

%compute 1/2*gradient

g = fjac' * ff;
