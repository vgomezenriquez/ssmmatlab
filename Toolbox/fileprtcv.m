function fileprtcv(hb, ttr, xx, tt, sconp)
%*************************************************************************
% This function prints the estimation results obtained with usa4vcvf.m
% Auxiliary function called in usa4vcv_d.m.
% Reference: ``Estimating Potential Output, Core Inflation
% and the NAIRU as Latent Variables'', by Rafael Domenech
% and Victor Gomez, Journal of Business and Economic Statistics (2006)
%
%   INPUTS:
%      hb : the beta estimate
%     ttr : t-values of the beta estimate
%      xx : estimated parameters
%      tt : t-values of the estimated parameters
%   sconp : standard deviation that has been concentrated out of the
%           likelihood
%
% Copyright (c) 21 July 2003 by Victor Gomez
% Ministerio de Hacienda, Direccion Gral. de Presupuestos,
% Subdireccion Gral. de Analisis y P.E.,
% Alberto Alcocer 2, 1-P, D-34, 28046, Madrid, SPAIN.
% E-mail: VGomez@sepg.minhap.es
%
% The author assumes no responsibility for errors or damage resulting from the use
% of this code. Usage of this code in applications and/or alterations of it should
% be referenced. This code may be redistributed if nothing has been added or
% removed and no money is charged. Positive or negative feedback would be appreciated.
%*************************************************************************


fid = fopen(fullfile('results', 'usa4vcv.out'), 'w');
fprintf(fid, '%s\n\n', 'Outliers:');
fprintf(fid, '%s%12.8f\n', '       O_1= ', hb(2));
fprintf(fid, '%s%6.2f%s\n', '              (', ttr(2), ')');
fprintf(fid, '%s%12.8f\n', '       O_2= ', hb(3));
fprintf(fid, '%s%6.2f%s\n', '              (', ttr(3), ')');
fprintf(fid, '%s%12.8f\n', '       O_3= ', hb(4));
fprintf(fid, '%s%6.2f%s\n', '              (', ttr(4), ')');
% fprintf(fid,'%s%12.8f\n','       O_4= ',hb(5));
% fprintf(fid,'%s%6.2f%s\n','              (',ttr(5),')');

fprintf(fid, '%s\n\n', 'Output:');
fprintf(fid, '%s%12.8f\n', '       g_y= ', hb(1));
fprintf(fid, '%s%6.2f%s\n', '              (', ttr(1), ')');
fprintf(fid, '%s%12.8f\n', 'sigma_{gw}= ', abs(xx(1)*sconp));
fprintf(fid, '%s%6.2f%s\n', '              (', tt(1), ')');
fprintf(fid, '%s%12.8f\n', '    theta1= ', xx(20));
fprintf(fid, '%s%6.2f%s\n', '              (', tt(20), ')');
fprintf(fid, '%s%12.8f\n', '    theta2= ', xx(21));
fprintf(fid, '%s%6.2f%s\n', '              (', tt(21), ')');
fprintf(fid, '%s%12.8f\n', 'sigma_{yw1}= ', abs(xx(5)*sconp));
fprintf(fid, '%s%6.2f%s\n', '              (', tt(5), ')');
fprintf(fid, '%s%12.8f\n', 'sigma_{yw2}= ', abs(xx(22)*sconp));
fprintf(fid, '%s%6.2f%s\n', '              (', tt(22), ')');

fprintf(fid, '%s\n\n', 'Unemployment:');
fprintf(fid, '%s%12.8f\n', 'sigma_{uw}= ', abs(xx(2)*sconp));
fprintf(fid, '%s%6.2f%s\n', '              (', tt(2), ')');
fprintf(fid, '%s%12.8f\n', '      phiu= ', xx(8));
fprintf(fid, '%s%6.2f%s\n', '              (', tt(8), ')');
fprintf(fid, '%s%12.8f\n', '      phi0= ', xx(14));
fprintf(fid, '%s%6.2f%s\n', '              (', tt(14), ')');
fprintf(fid, '%s%12.8f\n', '      phi1= ', xx(15));
fprintf(fid, '%s%6.2f%s\n', '              (', tt(15), ')');
fprintf(fid, '%s%12.8f\n', '      phi2= ', xx(16));
fprintf(fid, '%s%6.2f%s\n', '              (', tt(16), ')');
fprintf(fid, '%s%12.8f\n', 'sigma_{uv}= ', abs(xx(6)*sconp));
fprintf(fid, '%s%6.2f%s\n', '              (', tt(6), ')');

fprintf(fid, '%s\n\n', 'Investment:');
fprintf(fid, '%s%12.8f\n', 'sigma_{xw}= ', abs(xx(3)*sconp));
fprintf(fid, '%s%6.2f%s\n', '              (', tt(3), ')');
fprintf(fid, '%s%12.8f\n', '     betai= ', xx(9));
fprintf(fid, '%s%6.2f%s\n', '              (', tt(9), ')');
fprintf(fid, '%s%12.8f\n', '     beta0= ', xx(17));
fprintf(fid, '%s%6.2f%s\n', '              (', tt(17), ')');
fprintf(fid, '%s%12.8f\n', '     beta1= ', xx(18));
fprintf(fid, '%s%6.2f%s\n', '              (', tt(18), ')');
fprintf(fid, '%s%12.8f\n', 'sigma_{xv}= ', abs(xx(7)*sconp));
fprintf(fid, '%s%6.2f%s\n', '              (', tt(7), ')');

fprintf(fid, '%s\n\n', 'Inflation:');
fprintf(fid, '%s%12.8f\n', 'sigma_{pw1}=', abs(xx(4)*sconp));
fprintf(fid, '%s%6.2f%s\n', '              (', tt(4), ')');
fprintf(fid, '%s%12.8f\n', 'sigma_{pw2}=', abs(xx(23)*sconp));
fprintf(fid, '%s%6.2f%s\n', '              (', tt(23), ')');
fprintf(fid, '%s%12.8f\n', 'sigma_{pw3}=', abs(xx(24)*sconp));
fprintf(fid, '%s%6.2f%s\n', '              (', tt(24), ')');
fprintf(fid, '%s%12.8f\n', '        mu1=', xx(10));
fprintf(fid, '%s%6.2f%s\n', '              (', tt(10), ')');
fprintf(fid, '%s%12.8f\n', '        mu2=', xx(11));
fprintf(fid, '%s%6.2f%s\n', '              (', tt(11), ')');
fprintf(fid, '%s%12.8f\n', '        mu3=', xx(12));
fprintf(fid, '%s%6.2f%s\n', '              (', tt(12), ')');
fprintf(fid, '%s%12.8f\n', '        mu4=', xx(13));
fprintf(fid, '%s%6.2f%s\n', '              (', tt(13), ')');
fprintf(fid, '%s%12.8f\n', '       etay=', xx(19));
fprintf(fid, '%s%6.2f%s\n', '              (', tt(19), ')');
fprintf(fid, '%s%12.8f', 'sigma_{pv}= ', sconp);

fclose(fid);
