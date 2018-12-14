function tsplot(y, cstruc, varargin)
% PURPOSE: time-series plot with dates and labels
%---------------------------------------------------
% USAGE:     tsplot(y,cstruc,begp,endp,vnames,ydigit)
%        or: tsplot(y,cal_struc,vnames), which plots the entire series
%            with the default date format
%        or: tsplot(y,cal_struc,vnames,ydigit), which plots the entire series
%            with the date format ydigit
%        or: tsplot(y,cal_struc), entire series with no variable names
%        or: tsplot(y,cal_struc,[],ydigit), entire series with the date
%            format ydigit but no variable names
%
% where:  y      = matrix (or vector) of series to be plotted
%         cstruc = a structure returned by cal()
%         begp   = the beginning observation to plot (optional)
%         endp   = the ending observation to plot (optional)
%        vnames  = a string matrix of names for a legend (optional)
%                 e.g. vnames = ['y    ',
%                                'x1   ',  NOTE: fixed width
%                                'x2   ',        like all MATLAB
%                                'cterm'];       strings
%        ydigit  = a string specifying the date format (optional)
%                 e.g. ydigit = 'yyyy'
%                      ydigit = 'mmmyyyy'
%---------------------------------------------------
% e.g.   cstr = cal(1970,1,12);
%        tsplot(y,cstr); would plot all data
%
%    or: tsplot(y,cstr,ical(1981,1,cstr),ical(1981,12,cstr)),
%         which would plot data for 1981 only
%---------------------------------------------------
% SEE ALSO: tsprint
%---------------------------------------------------

% Original version of tsplot written by:
% James P. LeSage, Dept of Economics
% University of Toledo
% 2801 W. Bancroft St,
% Toledo, OH 43606
% jpl@jpl.econ.utoledo.edu
%
% Modified by Martyna Marczak, 20.09.2012
% Department of Economics (520G)
% University of Hohenheim
% Schloss, Museumsfluegel
% 70593 Stuttgart, Germany
% Phone: + 49 711 459 23823
% E-mail: marczak@uni-hohenheim.de
%**************************************************************************


nargs = length(varargin);

% Check of the arguments of tsplot corrected by Martyna Marczak,
% 20.09.2012
if nargin >= 2
    [nobs, nvar] = size(y);
    if ~isstruct(cstruc)
        error('tsplot: requires a cal() structure as input');
    end;
end;

freq = cstruc.freq;

if nargs == 0 % no user-supplied vnames or dates
    begp = 1;
    endp = nobs;
    if freq == 1
        ydigit = 'yyyy';
    elseif freq == 4
        ydigit = 'QQ-YY';
    elseif freq == 12
        ydigit = 'mmmyy';
    end
    vnames = [];
elseif nargs == 1
    if ischar(varargin{1}) || iscellstr(varargin{1}) % no dates but vnames
        begp = 1;
        endp = nobs;
        vnames = varargin{1};
        if freq == 1
            ydigit = 'yyyy';
        elseif freq == 4
            ydigit = 'QQ-YY';
        elseif freq == 12
            ydigit = 'mmmyy';
        end
    else
        error('The third argument must be a string/cell array if there are three arguments in tsplot');
    end
elseif nargs == 2; % we have an error
    % Why an error? begp and endp could be supplied, but vnames not
    if isscalar(varargin{1}) && isscalar(varargin{2})
        begp = varargin{1};
        endp = varargin{2};
        if freq == 1
            ydigit = 'yyyy';
        elseif freq == 4
            ydigit = 'QQ-YY';
        elseif freq == 12
            ydigit = 'mmmyy';
        end
        vnames = [];
    elseif (ischar(varargin{1}) || iscellstr(varargin{1}) || ...
            isempty(varargin{1})) && ischar(varargin{2})
        begp = 1;
        endp = nobs;
        vnames = varargin{1};
        ydigit = varargin{2};
    else
        error('The third and the fourth argument must be either both number or a string/cell array/empty and a string, respectively, if there are four arguments in tsplot');
    end
elseif nargs == 3
    if isscalar(varargin{1}) && isscalar(varargin{2})
        if ischar(varargin{3}) || iscellstr(varargin{3})
            begp = varargin{1};
            endp = varargin{2};
            vnames = varargin{3};
            if freq == 1
                ydigit = 'yyyy';
            elseif freq == 4
                ydigit = 'QQ-YY';
            elseif freq == 12
                ydigit = 'mmmyy';
            end
        else
            error('The fifth argument must be a string/cell array if there are five arguments in tsplot');
        end
    else
        error('The third and the fourth argument must be numbers if there are five arguments in tsplot')
    end
elseif nargs == 4
    if isscalar(varargin{1}) && isscalar(varargin{2})
        if (ischar(varargin{3}) || iscellstr(varargin{3}) || isempty(varargin{3}))
            if ischar(varargin{4})
                begp = varargin{1};
                endp = varargin{2};
                vnames = varargin{3};
                ydigit = varargin{4};
            else
                error('The sixth argument must be a string if there are six arguments in tsplot');
            end
        else
            error('The fifth argument must be a string/cell array or empty if there are six arguments in tsplot')
        end
    else
        error('The third and the fourth argument must be numbers if there are five arguments in tsplot');
    end
end

if isempty(vnames) % no variable names supplied, make some up
    for i = 1:nvar
        if i < 10
            snames = 'series  ';
            name = [snames, num2str(i)];
            vnames = [vnames, name];
        else
            snames = 'series ';
            name = [snames, num2str(i)];
            vnames = [vnames, name];
        end
    end
else
    [vsize, nsize] = size(vnames); % error check vnames argument
    if vsize ~= nvar
        error('Wrong # vnames in tsplot');
    end
end


fsize = 9; % font size
[nobs, nvar] = size(y(begp:endp, :)); % find nobs, nvar


if nobs <= 120; % provide a grid for small samples
    grid = 'on';
else
    grid = 'off';
end;

%******************************************************
% Extension of the plots by Martyna Marczak, 20.09.2012
switch freq;
    case 1, % case of annual series
        out = cal(cstruc.beg_yr, cstruc.beg_per, cstruc.freq, begp);
        beg_yr = out.year;
        yr = beg_yr:beg_yr + nobs - 1;
        yrs = yr';
        % Graphics handle added
        g = plot(datenum(yrs, 1, 1), y(begp:endp, :));
        legend(g, vnames, 'Location', 'NorthWest'); % least conflict with the data plot
        % Add a line at y=0 if some of the y values lie below 0
        for i = nvar
            if any(y(begp:endp, i) < 0)
                hold on
                xlim = get(gca, 'XLim');
                plot(xlim, [0, 0], ':', 'Color', [0.4, 0.4, 0.4], ...
                    'LineWidth', 0.05)
                hold off
                break
            end
        end
    case 4, % case of quarterly series
        yrs = zeros(nobs, 1);
        qtr = zeros(nobs, 1);
        out = cal(cstruc.beg_yr, cstruc.beg_per, cstruc.freq, begp);
        beg_yr = out.year;
        %beg_qtr = out.period;
        % BUG fix suggested by Stephen Burke
        % PhD Student
        % Faculty of Commerce and Bus. Admin, Dept. of Finance
        % University of British Columbia
        if out.period == 1
            %beg_qtr = 1;
            beg_qtr = 3;
        elseif out.period == 2
            %beg_qtr = 4;
            beg_qtr = 6;
        elseif out.period == 3
            %beg_qtr = 7;
            beg_qtr = 9;
        else
            %beg_qtr = 10;
            beg_qtr = 12;
        end
        for i = 1:nobs;
            yrs(i, 1) = beg_yr;
            qtr(i, 1) = beg_qtr;
            beg_qtr = beg_qtr + 3;
            if beg_qtr > 12
                beg_yr = beg_yr + 1;
                beg_qtr = 1;
            end
        end
        % Graphics handle added
        g = plot(datenum(yrs, qtr, 1), y(begp:endp, :));
        legend(g, vnames, 'Location', 'NorthWest'); % least conflict with the data plot
        % Add a line at y=0 if some of the y values lie below 0
        for i = nvar
            if any(y(begp:endp, i) < 0)
                hold on
                xlim = get(gca, 'XLim');
                plot(xlim, [0, 0], ':', 'Color', [0.4, 0.4, 0.4], ...
                    'LineWidth', 0.05)
                hold off
                break
            end
        end
    case 12, % case of monthly series
        yrs = zeros(nobs, 1);
        mth = zeros(nobs, 1);
        out = cal(cstruc.beg_yr, cstruc.beg_per, cstruc.freq, begp);
        %          if strcmpi(ydigit,'yyyy')
        %             beg_yr = out.year+1;
        %          else
        %            beg_yr = out.year;
        %          end
        beg_yr = out.year;
        beg_mth = out.period;
        for i = 1:nobs;
            yrs(i, 1) = beg_yr;
            mth(i, 1) = beg_mth;
            beg_mth = beg_mth + 1;
            if beg_mth > 12
                beg_yr = beg_yr + 1;
                beg_mth = 1;
            end
        end
        g = plot(datenum(yrs, mth, 1), y(begp:endp, :));
        legend(g, vnames, 'Location', 'NorthWest'); % least conflict with the data plot
        % Add a line at y=0 if some of the y values lie below 0
        for i = nvar
            if any(y(begp:endp, i) < 0)
                hold on
                xlim = get(gca, 'XLim');
                plot(xlim, [0, 0], ':', 'Color', [0.4, 0.4, 0.4], ...
                    'LineWidth', 0.05)
                hold off
                break
            end
        end
    otherwise % how did we get here?
        disp('frequency unknown to tsplot');
end

%*******************************************************
% Extension by Martyna Marczak, 20.09.2012
if yrs(end) - yrs(1) < 70
    yrsb = yrs(1) / 10;
    yrse = yrs(end) / 10;
    diffs = (yrsb - floor(yrsb)) * 10;
    diffe = (ceil(yrse) - yrse) * 10;
    if yrs(end) - yrs(1) < 10
        % Make the ticks each year instead of 10 years
        nb = 1;
        yrsst = yrs(1);
        yrsend = yrs(end);
    elseif (yrs(end) - yrs(1) >= 10) && (yrs(end) - yrs(1) < 30)
        % Make the ticks each two years instead of 10 years
        nb = 2;
        if diffs < 2
            yrsst = floor(yrsb) * 10;
        else
            yrsst = floor(yrsb) * 10 + 2;
        end
        if diffe < 2
            yrsend = ceil(yrse) * 10;
        else
            yrsend = ceil(yrse) * 10 - 2;
        end
    else
        % Make the ticks each 5 years instead of 10 years
        nb = 5;
        if diffs < 5
            yrsst = floor(yrsb) * 10; % start year's last number is 0
        else
            yrsst = floor(yrsb) * 10 + 5; % start year's last number is 5
        end
        if diffe < 5
            yrsend = ceil(yrse) * 10; % end years's last number is 0
        else
            yrsend = ceil(yrse) * 10 - 5; % end years's last number is 5
        end
    end
    nticks = ((yrsend - yrsst) / nb) + 1;
    xticks = zeros(1, nticks);
    nb1 = nb;
    for i = 1:nticks
        xticks(i) = datenum(yrsst+nb1, 1, 1);
        nb1 = nb1 + nb;
    end
    set(gca, 'Xtick', xticks)
end
%*********************************************************
% Axes settings

set(gca, 'fontsize', fsize);
set(gca, 'tickdir', 'in');
datetick('x', ydigit, 'keepticks');
set(gca, 'GridLineStyle', ':');
%set(gca,'xgrid',grid);       %comment this line if grid is not desired
set(gca, 'xcolor', 'k');

% Set y axis limits

ymin = min(min(y));
% Count the number of digits
if ymin ~= 0
    dig1 = real(floor(log10(ymin)+1));
else
    dig1 = 1;
end
if dig1 <= 0
    f1 = 0;
    while dig1 <= 0
        f1 = f1 + 1;
        ff1 = 10^f1;
        ymin1 = ymin * ff1;
        dig1 = real(floor(log10(ymin1)+1));
    end
    f1 = 1 / f1;
    fact1 = 1 / (10^f1);
else
    fact1 = 1 / (10^dig1);
end
ylim1 = ymin - fact1 * abs(ymin);


ymax = max(max(y));
if ymax ~= 0
    dig2 = real(floor(log10(ymax)+1));
else
    dig2 = 1;
end
if dig2 <= 0
    f2 = 0;
    while dig2 <= 0
        f2 = f2 + 1;
        ff2 = 10^f2;
        ymax1 = ymax * ff2;
        dig2 = real(floor(log10(ymax1)+1));
    end
    f2 = 1 / f2;
    fact2 = 1 / (10^f2);
else
    fact2 = 1 / (10^dig2);
end
ylim2 = ymax + fact2 * abs(ymax);
set(gca, 'YLim', [ylim1, ylim2]);


% Set x axis limits
if freq == 1
    xlim1 = datenum(yrs(1), 1, 1);
    xlim2 = datenum(yrs(end), 1, 1);
elseif freq == 4
    xlim1 = datenum(yrs(1), qtr(1), 1);
    xlim2 = datenum(yrs(end), qtr(end), 1);
else
    xlim1 = datenum(yrs(1), mth(1), 1);
    xlim2 = datenum(yrs(end), mth(end), 1);
end

set(gca, 'XLim', [xlim1, xlim2]);
