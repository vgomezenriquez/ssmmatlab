function plotres(y, Y, g, yor, datei, cw, fname, gft, nrout, Youtg, nreg, Yrg, infr, s, lam)
%**************************************************************************
% This function plots original series, residuals, outliers, regression
% variables, residual histogram, and correlograms of residuals and squared
% residuals.
%
%     INPUTS:
%         y : data vector
%         Y : matrix with regression variables
%         g : vector with regression coeffients
%       yor : original time series
%     datei : calendar structure
%        cw : critical value of the standard normal distribution to compute
%             confidence bounds
%     fname : series label appearing in the legend
%       gft : =1 : show and do not save figures
%             =2 : save and do not show figures
%             =3 : save and show figures
%     nrout = 0 : do not produce graph of outlier effects
%           > 0 : produce graph for each outlier effect
%     Youtg : matrix with outlier effects
%      nreg = 0 : do not produce graph for the effects of regression
%                 variables other than outliers
%           > 0 : produce graph for the effects of each regression variable
%                 other than outlier
%       Yrg : matrix with effects of regression variables other than
%             outliers
%      infr : residual structure (output of rescomp)
%         s : frequency of the data
%       lam = 0 : compute logs of the original series
%           = 1 : do not compute logs
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
%**************************************************************************

if ~isstruct(infr)
    error('plotres: requires a residual structure');
end;
if ~isstruct(datei) && ~isempty(y)
    error('plotres: requires a calendar structure');
end;

if isfield(infr, 'ycii')
    ycii = infr.ycii;
else
    ycii = [];
end

e = infr.e;
ne = infr.ne;
stde = infr.stde;
r = infr.r;
pc = infr.pc;
sea = infr.sea;
sep = infr.sep;
no = infr.no;
xo = infr.xo;
rs = infr.rs;
pcs = infr.pcs;
seas = infr.seas;

if (gft == 3) || (gft == 1)
    gflag = 1;
else
    gflag = 0;
end
pathc = pwd; %current directory
if ~isempty(y)
    ny = length(y);
    if gft == 2
        f = figure('visible', 'off');
    else
        f = figure('visible', 'on');
    end
    if (s == 12) || (s == 4)
        tsplot(yor, datei, fname)
    else
        plot(yor);
        legend(fname);
    end
    if (gft > 1)
        saveas(f, [pathc, filesep, 'graphs', filesep, 'Original'], 'pdf')
    end
    if (gflag == 1)
        disp('strike any key when ready')
        pause
        close all
    end
    flagycii = 0;
    if ~isempty(ycii)
        if gft == 2
            f = figure('visible', 'off');
        else
            f = figure('visible', 'on');
        end
        if lam == 0
            vnames = char('Original Series in logs','Series corrected by filtered inputs only');
        else
            vnames = char('Original Series','Series corrected by filtered inputs only');
        end
        if (s == 12) || (s == 4)
            tsplot([y, ycii], datei, vnames);
        else
            t = 1:ny;
            plot(t, y, t, ycii);
            legend(vnames);
        end
        if (gft > 1)
            saveas(f, [pathc, filesep, 'graphs', filesep, 'Originalcorrfi'], 'pdf')
        end
        if (gflag == 1)
            disp('strike any key when ready')
            pause
            close all
        end
        flagycii = 1;
    end
    [ng, mg] = size(g);
    if ng > 0
        if gft == 2
            f = figure('visible', 'off');
        else
            f = figure('visible', 'on');
        end
        if flagycii == 0
            ycor = y - Y(1:ny, :) * g;
        else
            ycor = ycii - Y(1:ny, :) * g;
        end
        if lam == 0
            vnames = char('Original Series in logs','Corrected Series');
        else
            vnames = char('Original Series','Corrected Series');
        end
        if (s == 12) || (s == 4)
            tsplot([y, ycor], datei, vnames);
        else
            t = 1:ny;
            plot(t, y, t, ycor);
            legend(vnames);
        end
        if (gft > 1)
            saveas(f, [pathc, filesep, 'graphs', filesep, 'Originalcorr'], 'pdf')
        end
        if (gflag == 1)
            disp('strike any key when ready')
            pause
            close all
        end
    end
end

if nreg > 0
    if gft == 2
        f = figure('visible', 'off');
    else
        f = figure('visible', 'on');
    end
    [nYr, mYr] = size(Yrg);
    t = 1:nYr;
    %  for i=1:nreg
    %   plot(t,Yrg(:,i)); hold on
    %  end
    Yrgt = zeros(nYr, 1);
    for i = 1:nYr
        Yrgt(i) = sum(Yrg(i, 1:nreg));
    end
    plot(t, Yrgt);
    title('Sum of regression effects other than outliers'); %hold off;
    if (gft > 1)
        saveas(f, [pathc, filesep, 'graphs', filesep, 'Sumofregef'], 'pdf')
    end
    if (gflag == 1)
        disp('strike any key when ready')
        pause
        close all
    end
end

if nrout > 0
    if gft == 2
        f = figure('visible', 'off');
    else
        f = figure('visible', 'on');
    end
    [nYo, mYo] = size(Youtg);
    t = 1:nYo;
    for i = 1:nrout
        plot(t, Youtg(:, i));
        hold on
    end
    title('Outlier effects');
    hold off;
    if (gft > 1)
        saveas(f, [pathc, filesep, 'graphs', filesep, 'Outliereff'], 'pdf')
    end
    if (gflag == 1)
        disp('strike any key when ready')
        pause
        close all
    end
end

if gft == 2
    f = figure('visible', 'off');
else
    f = figure('visible', 'on');
end
t = 1:ne;
bd = ones(ne, 1) * (cw * stde); %confidence bands
plot(t, e, t, zeros(ne, 1), t, bd, 'r', t, -bd, 'r') %plot residuals
legend('Residuals')
if (gft > 1)
    saveas(f, [pathc, filesep, 'graphs', filesep, 'Residuals'], 'pdf')
end
if (gflag == 1)
    disp('strike any key when ready')
    pause
    close all
end

if gft == 2
    f = figure('visible', 'off');
else
    f = figure('visible', 'on');
end
bar(xo, no);
title('Residual histogram')
xlabel('Standard deviation intervals')
if (gft > 1)
    saveas(f, [pathc, filesep, 'graphs', filesep, 'Residualhist'], 'pdf')
end
if (gflag == 1)
    disp('strike any key when ready')
    pause
    close all
end


if gft == 2
    f = figure('visible', 'off');
else
    f = figure('visible', 'on');
end
fname = 'residuals';
rpplot(r, pc, sea, sep, cw, fname); %plot residual autocorrelations
if (gft > 1)
    saveas(f, [pathc, filesep, 'graphs', filesep, 'Resoutcor'], 'pdf')
end
if (gflag == 1)
    disp('strike any key when ready')
    pause
    close all
end


if gft == 2
    f = figure('visible', 'off');
else
    f = figure('visible', 'on');
end
fname = 'squared residuals';
rpplot(rs, pcs, seas, sep, cw, fname); %plot autocorrelations
if (gft > 1)
    saveas(f, [pathc, filesep, 'graphs', filesep, 'SqResoutcor'], 'pdf')
end
if (gflag == 1)
    disp('strike any key when ready')
    pause
    close all
end
