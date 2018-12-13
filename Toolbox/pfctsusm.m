function pfctsusm(out)
%*************************************************************************
% This function plots forecasts
%
%      INPUT:
%       out : structure with the following fields:
%      .yor : original time series
%        .y : time series used in computation of forecasts
%       .ny : length of y
%      .pry : forecasts of y
%     .spry : standard errors of pry
%     .opry : forecasts of y in the original scale (if lam = 0)
%    .ospry : standard errors of opry
%      .npr : number of forecasts
%       .cw : critical value of the standard normal distribution used in
%             computation of confidence bounds
%    .tname : name of the time series
%        .s : frequency of the data
%      .lam = 0 : compute logs of yor
%           = 1 : do not compute logs of yor
%      .gft : =1 : show and do not save figures
%             =2 : save and do not show figures
%             =3 : save and show figures
%
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

if ~isfield(out, 'gft')
    gft = 1;
else
    gft = out.gft;
end
if (gft == 3) || (gft == 1)
    gflag = 1;
else
    gflag = 0;
end
pathc = pwd; %current directory
pry = out.pry;
spry = out.spry;
opry = out.opry;
ospry = out.ospry;
y = out.y;
yor = out.yor;
ny = out.ny;
npr = out.npr;
cw = out.cw;
tname = out.tname;
lam = out.lam;
s = out.s;
if gft == 2
    f = figure('visible', 'off');
else
    f = figure('visible', 'on');
end
tt = ny - npr + 1:ny;
pry = pry';
spry = spry';
opry = opry';
ospry = ospry';
y1 = [y(tt); pry + cw * spry];
y2 = [y(tt); pry];
y3 = [y(tt); pry - cw * spry];
t = -npr + 1:npr;
vnames = ['Upper 95% band '; ...
    'Forecast       '; ...
    'Lower 95% band '];
plot(t, y1, 'r-.', t, y2, '-', t, y3, '--');
legend(vnames);
set(gca, 'tickdir', 'in');
set(gca, 'xcolor', 'k');
set(gca, 'GridLineStyle', ':');
set(gca, 'Xgrid', 'on');
axis tight;

if lam == 0
    title(['Forecasts of series ', tname, ' (in logs)'])
    if (gft > 1)
        saveas(f, [pathc, filesep, 'graphs', filesep, 'Forecastl'], 'pdf')
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
    y1 = [yor(tt); opry + cw * sqrt(ospry)];
    y2 = [yor(tt); opry];
    y3 = [yor(tt); opry - cw * sqrt(ospry)];
    t = -npr + 1:npr;
    vnames = ['Upper 95% band '; ...
        'Forecast       '; ...
        'Lower 95% band '];
    plot(t, y1, 'r-.', t, y2, '-', t, y3, '--');
    legend(vnames);
    set(gca, 'tickdir', 'in');
    set(gca, 'xcolor', 'k');
    set(gca, 'GridLineStyle', ':');
    set(gca, 'Xgrid', 'on');
    axis tight;
    title(['Forecasts of series ', tname])
    if (gft > 1)
        saveas(f, [pathc, filesep, 'graphs', filesep, 'Forecast'], 'pdf')
    end
    if (gflag == 1)
        disp('strike any key when ready')
        pause
        close all
    end
else
    title(['Forecasts of series  ', tname])
    if (gft > 1)
        saveas(f, [pathc, filesep, 'graphs', filesep, 'Forecast'], 'pdf')
    end
    if (gflag == 1)
        disp('strike any key when ready')
        pause
        close all
    end
end
if (s > 1)
    %plot forecasts in rates of grothw
    if gft == 2
        f = figure('visible', 'off');
    else
        f = figure('visible', 'on');
    end
    yt = tasa([yor; opry], s) * 100;
    nyt = length(yt);
    tt = nyt - 2 * npr + 1:nyt;
    t = -npr + 1:npr;
    plot(t, yt(tt), 'r-');
    set(gca, 'tickdir', 'in');
    set(gca, 'xcolor', 'k');
    set(gca, 'GridLineStyle', ':');
    set(gca, 'Xgrid', 'on');
    axis tight;
    title(['Forecasts of original series ', tname, ' (rates of growth in percentage)'])
    if (gft > 1)
        saveas(f, [pathc, filesep, 'graphs', filesep, 'Forecastrg'], 'pdf')
    end
    if (gflag == 1)
        disp('strike any key when ready')
        pause
        close all
    end
end
