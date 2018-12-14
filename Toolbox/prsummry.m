function prsummry(ii, ny, nreg0, fid, fname, parm, iout)
%
% this function prints the automatic model identification summary for a
% list of ARIMA or transfer function models.
%
%     INPUTS:
%     ii     : an integer, corresponding to the series currently handled.
%     ny     : the series length
%     nreg0  : the number of original regression variables
%     fid    : the number of the output file
%     fname  : a string containing the series name
% parm: astructure containing model infomation, where
% .s:  seasonality
% .p:  AR order
% .ps: order of the AR of order s
% .q:  order of the regular MA
% .qs: order of the MA of order s (1 at most)
% .dr: order of regular differencing
% .ds: order of differencing of order s
% .lam: = 0, logs are taken, = 1, no logs
%.flagm: = 0, no mean, = 1, mean in the model
%.trad : the estimated number of TD variables
%.leap : = 0, no leap year effect, = 1, leap year effect in the model
%.east : = 0, no Easter effect, = 1, Easter effect in the model
%.dur  : duration of the Easter effect
% iout:   a structure containing information for outlier detection, where
% .C:     critical value for outlier detection
% .delta: the value for delta in TC outliers
% .mthd: method to compute ARMA parameter estimates (0 Hannan Rissanen, 1 Max. Lik.)
% .schr: =0 outliers of type AO and TC are considered (default)
%        =1 outliers of type AO, TC and LS are considered
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

p = parm.p;
dr = parm.dr;
q = parm.q;
ps = parm.ps;
ds = parm.ds;
qs = parm.qs;
s = parm.s;
lam = parm.lam;
flagm = parm.flagm;
trad = parm.trad;
leap = parm.leap;
east = parm.east;
dur = parm.dur;
nrout = iout.nrout;
C = iout.C;
omet = iout.omet;
nind = iout.nind;
tip = iout.tip;
fprintf(fid, '%s %1i %1s %1i %1s %1i %s', [fname, ' ('], p, ',', dr, ',', q, ')');
if s > 1
    fprintf(fid, '%1s %1i %1s %1i %1s %1i %s %2i', '(', ps, ',', ds, ',', qs, ')_', s);
end
fprintf(fid, ' lam= %1i', lam);
fprintf(fid, ' mean= %1i', flagm);
fprintf(fid, ' nreg0= %2i', nreg0);
fprintf(fid, ' trade= %1i', trad);
fprintf(fid, ' leapy= %1i', leap);
fprintf(fid, '\n');
fprintf(fid, '   %5i', ii);
fprintf(fid, '  ny= %i', ny);
fprintf(fid, [' east(', num2str(dur), ')= %1i'], east);
fprintf(fid, '  nrout= %i', nrout);
fprintf(fid, '  C= %3.1f', C);
if omet == 1, met = 'EML';
else, met = 'HR ';
end
fprintf(fid, '  met= %s', met);
if nrout > 0
    count = nrout;
    counti = 0;
    outcol = 9;
    while count > 0
        fprintf(fid, '\n');
        if counti == 0, fprintf(fid, '          out= '); ...
        else, fprintf(fid, '             ');
        end
        for i = 1:min(outcol, count)
            fprintf(fid, '%4i  %2s', nind(i+counti), tip(i+counti, 1:2));
        end
        count = count - outcol;
        counti = counti + i - 1;
    end
end
fprintf(fid, '\n');
