function trtout(fid,iout,ser)
%
% This function prints the results of outlier detection
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

C=iout.C; omet=iout.omet; nind=iout.nind; tip=iout.tip; ornames=iout.ornames;
nrout=iout.nrout; freq=ser.freq; bg_year=ser.bg_year; bg_per=ser.bg_per;
if omet == 1, met='Exact max. likelihood'; else, met='Hannan-Rissanen'; end
fprintf(fid,'%23s %3.1f %s\n','Outliers detected (C = ',C,[', Method is ' met '):']); %print header
fprintf(fid,'%6s  %11s  %4s  %4s  %6s\n','Order','Obs. number',...
 'Type','Year','Period');
for i=1:nrout
 st=cal(bg_year,bg_per,freq,nind(i)); 
 fprintf(fid,'%6s  %11.f  %4s  %4.f  %6.f\n',ornames(i,:),...
  nind(i),tip(i,1:2),st.year,st.period);
end
