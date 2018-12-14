function nberrecplot(recdates, y, color)
%**************************************************************************
%    Function nberrecplot plots areas specifying recession dates
%
%       INPUTS:
%      REQUIRED
%   recdates  : (mrec x 4) matrix with the peak and trough dates;
%                mrec is the number of recessions;
%               Columns:
%               1. year of each peak
%               2. month of year (from 1.) of each peak
%               3. year of each trough
%               4. month of year (from 3.) of each trough
%          y  : vector with at least two elements;
%               y can be a series or a vector with two elements,
%               in each case the minimum and the maximum value will be
%               used for specifying the height of the rectangular area;
%
%     OPTIONAL(the order does not matter)
%      color  : color specification, e.g.
%               color = 'red'
%               color = [0.6,0.5,0.9];
%               default is [0.8,0.8,0.8] (gray)
%**************************************************************************
% Written by: Martyna Marczak, 20.09.2012
% Department of Economics (520G)
% University of Hohenheim
% Schloss, Museumsfluegel
% 70593 Stuttgart, Germany
% Phone: + 49 711 459 23823
% E-mail: marczak@uni-hohenheim.de
%**************************************************************************

% Check the number and order of the arguments and set defaults

if nargin < 2
    error('There must be at least two inputs to nberrecplot');
end

[mrec, nrec] = size(recdates);

if nrec ~= 4
    error('recdates must have four columns');
else
    if iscell(recdates)
        error('recdates must be a matrix');
    end
end

if ~isvector(y)
    error('y must be a vector')
else
    if length(y) < 2
        error('y must have at least two elements');
    end
end

if nargin == 2
    color = [0.8, 0.8, 0.8];
end

if nargin == 3
    if isvector(color) && ((~ischar(color)) && ~iscellstr(color))
        if length(color) ~= 3
            error('color must have three elements if it is a vector');
        else
            for i = 1:3
                if color(i) < 0 || color(i) > 1
                    error('Each element of color should take a value [0,1]');
                end
                break
            end
        end
    elseif ischar(color)
        co = {'y', 'm', 'c', 'r', 'g', 'b', 'w', 'k'};
        col = {'yellow', 'magenta', 'cyan', 'red', 'green', 'blue', 'white', 'black'};
        if (~any(strcmpi(color, co))) && (~any(strcmpi(color, col)))
            error('color must be an admissible string');
        end
    else
        error('color must be either a vector or a string');
    end
end


ymin = min(y);
ylmin = ymin;
%ylmin = ymin - 0.09*abs(ymin);
ymax = max(y);
ylmax = ymax;
%ylmax = ymax + 0.09*abs(ymax);


% Convert the date to serial date number
pdatenum = zeros(mrec, 1);
tdatenum = zeros(mrec, 1);

for i = 1:mrec
    pdate = recdates(i, :); % date of the peak
    pdatenum(i) = datenum(pdate(1), pdate(2), 1);
    tdate = recdates(i, 3:4); % date of the trough
    tdatenum(i) = datenum(tdate(1), tdate(2), 1);
end

for i = 1:mrec
    h = fill([pdatenum(i), tdatenum(i), tdatenum(i), pdatenum(i)], ...
        [ylmin, ylmin, ylmax, ylmax], ...
        color, 'LineStyle', 'none');
    hold on
end
