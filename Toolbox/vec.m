function v = vec(x)
% PURPOSE: creates a column vector by stacking columns of x
%----------------------------------------------------------
% USAGE:  v = vec(x)
% where:  x = an input matrix
%---------------------------------------------------------
% RETURNS:
%         v = output vector containing stacked columns of x
%----------------------------------------------------------

% Written by KH (Kurt.Hornik@tuwien.ac.at) on 1995/05/08
% Copyright Dept of Probability Theory and Statistics TU Wien

% Modified by J.P. LeSage

  if (nargin ~= 1)
  error('Wrong # of arguments to vec'); 
  end 
  
  v = x(:);
  
  
