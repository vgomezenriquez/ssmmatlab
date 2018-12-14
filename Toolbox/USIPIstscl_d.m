%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% US Industrial Production Index, 1946-I, 2011-III
% A univariate structural model is specified and estimated. It includes a
% cyclical component and an outlier (LS). The outlier is assigned to the
% trend.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;

%model estimation
%the outlier is incorporated to the trend through the field Ycomp of
%structure ser in the specification file usmUSIPI
out = usmestos('usmUSIPI');
