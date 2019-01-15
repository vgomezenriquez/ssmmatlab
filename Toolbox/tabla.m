function tabla(y, datei, info)
%
% This function produces a table for a time series
%
% Input arguments:
%               y: the time series
%           datei: a structure containing the initial date
%            info: a structure containing printing options,
%  where
%            .fid: the device on which the table will be written
%          .fh: flag for header and years
%             .wd: format width
%             .nd: number of decimal points
%          .scale: =1 scale data if necessary
%                  =0 do not scale data
%
%    Note: this routine can be used to print data without a header and years
%          in several columns using datei for the number of columns.
% Example: datei=cal(1900,1,6) and info.fh=0 will print only the data with six
%          columns. The year 1900 is immaterial here.
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
%


if ~isstruct(datei)
    error('tabla: requires a calendar structure');
end;
% setup defaults
fid = 1;
fh = 1;
wd = 10;
nd = 4;
scale = 0;
if ~isempty(info)
    if ~isstruct(info)
        error('tabla: options should be in a structure variable');
    end;
    % parse options
    fields = fieldnames(info);
    nf = length(fields);
    for i = 1:nf
        if strcmp(fields{i}, 'fid')
            fid = info.fid;
        elseif strcmp(fields{i}, 'fh')
            fh = info.fh;
        elseif strcmp(fields{i}, 'wd')
            wd = info.wd;
        elseif strcmp(fields{i}, 'nd')
            nd = info.nd;
        elseif strcmp(fields{i}, 'scale')
            scale = info.scale;
        end;
    end;
else
    % rely on default options
end;

fsmall = 0;
fbig = 0;
if scale == 1
    eps = 1e-10;
    ys = y + ones(size(y)) * eps; %scale data if necessary
    minly = min(log10(abs(ys)));
    if minly < 0, minly = -minly;
    else minly = 0;
    end
    if minly > 1
        small = min(9, floor(minly));
        fsmall = 1;
        y = y * 10^small;
    end
    maxly = max(log10(abs(ys)));
    if maxly > 5
        big = min(9, floor(maxly)) - 3;
        fbig = 1;
        y = y * 10^(-big);
    end
end

datef = cal(datei, length(y)); %creates a structure containing the starting and final dates
fy = datef.beg_yr; %first year
ny = datef.year - fy + 1; %number of years
ns = datef.freq; %number of seasons
npi = datei.beg_per;
npf = datef.period; %initial and final periods
fwd = num2str(wd);
fnd = num2str(nd);
fmd = ['%', fwd, '.', fnd, 'f']; %format for the data
if fh == 1
    fmi = '%5i '; %format for the years
elseif fh == 0
    fmi = '%1s';
end
%write header
if (ns == 12)
    yearc = ['Jan'; 'Feb'; 'Mar'; 'Apr'; 'May'; 'Jun'; 'Jul'; 'Aug'; 'Sep'; 'Oct'; 'Nov'; 'Dec'];
elseif (ns > 1) & (ns <= 6)
    yearc = ['  I'; ' II'; 'III'; ' IV'; '  V'; ' VI'];
end
if (ns > 1) & ((ns == 12) | (ns == 6) | (ns == 5) | (ns == 4) | (ns == 3) | ...
        (ns == 2)) & (fh == 1) fmby = ['%', fwd, 's'];
    fprintf(fid, '%5s ', 'Year');
    for i = 1:ns, fprintf(fid, fmby, strjust(yearc(i, :), 'right'));
    end
    fprintf(fid, '\n');
end
%first row
t = 1:ns - (npi - 1); %observation indices
if fh == 1, fprintf(fid, fmi, fy);
elseif fh == 0, fprintf(fid, fmi, ' ');
end
if npi > 1
    nb = num2str(wd*(npi - 1));
    fmb = ['%', nb, 's']; %format for the first blanks
    fprintf(fid, fmb, ' ');
    fprintf(fid, fmd, y(t));
    fprintf(fid, '\n');
else
    fprintf(fid, fmd, y(t));
    fprintf(fid, '\n');
end
%intermediate rows
if ny > 1
    for i = 2:ny - 1
        fy = fy + 1;
        t = (i - 1) * ns - (npi - 2):i * ns - (npi - 1);
        if fh == 1, fprintf(fid, fmi, fy);
        elseif fh == 0, fprintf(fid, fmi, ' ');
        end
        fprintf(fid, fmd, y(t));
        fprintf(fid, '\n');
    end
end
%last row
if ny > 2
    fy = fy + 1;
    t = (ny - 1) * ns - (npi - 2):(ny - 1) * ns - npi + npf + 1;
    if fh == 1, fprintf(fid, fmi, fy);
    elseif fh == 0, fprintf(fid, fmi, ' ');
    end
    fprintf(fid, fmd, y(t));
    fprintf(fid, '%s\n', '');
end
%write note on scaling if necessary
if fsmall == 1, fprintf(fid, ' _________\n');
    fprintf(fid, ' Values should be multiplied by 1e-%i', small);
end
if fbig == 1, fprintf(fid, ' _________\n');
    fprintf(fid, ' Values should be multiplied by 1e+%i', big);
end
