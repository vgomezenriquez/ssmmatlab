function Y = trade(Iy,Im,N,Itrad,Ntrad,Yh,Mq)

% this function generates the trading day variables and the leap year variable.  
% It works until 2100.
%
% input variables       Iy      : the initial year
%                       Im      : the initial period
%                       N       : the length of the desired vector
%                       Itrad   : the number of trading day variables
%                                 =1 one trading day variable (sundays and
%                                 saturdays are the holidays)
%                                 =2 like 1, but the leap year variable is 
%                                    added
%                                 =6 six trading day variables are considered
%                                 =7 like 6, but the leap year variable is 
%                                    added
%                       Ntrad   : a flag for additional holydays (=1 yes, =0 no)
%                       Yh      : the variable where the additional holidays are 
%                                 stored
%                       Mq      : the series frequency (=12 for monthly, 
%                                                       =4 for quarterly)
%
% output variables      Y       : N x Itrad array containing the trading day 
%                                 variables%
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


Y=zeros(N,Itrad);
% quarterly series
      if (Mq == 4), in=N; Im=(Im-1)*3+1; N=N*3; end
      dy=floor((Im-1)/12); mm=Im-12*dy; yy=Iy+dy; ly=mod(yy,4);
% day of the week of the first day of the year
      y=yy-1901; l=floor(y/4);  d=y+l+2; dd=mod(d,7)+1;
% day of the week of the first day of month Im in year Iy
      if (mm == 5) 
       dd=dd+1;
      elseif (mm == 8)
       dd=dd+2;
      elseif (mm == 2) | (mm == 3) | (mm == 11)
       dd=dd+3;
      elseif (mm == 6)
       dd=dd+4;
      elseif (mm == 9) | (mm == 12) 
       dd=dd+5;
      elseif (mm == 4) | (mm == 7) 
       dd=dd+6;
      end
      if ((ly == 0) & (mm >= 3)), dd=dd+1; end
      if (dd > 7), dd=dd-7; end
% store in z0 the distribution of the days of the week corresponding to month i in the form
% sat fri thu wed tue mon sun
      z0=zeros(8,1);
      for i=1:N
       for j=1:7
        z0(j)=4.D0;
       end
       flag=0;
       if ( (mm == 1) | (mm == 3) | (mm == 5) | (mm == 7) | (mm == 8) | (mm == 10) | (mm == 12))
        dif=3; z0(8)=31.D0;
       elseif ((mm == 4) | (mm == 6) | (mm == 9) | (mm == 11)) 
        dif=2; z0(8)=30.D0;
       elseif (mm == 2) 
        dif=0; z0(8)=28.D0;
        if (ly ~=0), flag=1; else, dif=1; z0(8)=29.D0; end
       end
       if flag == 0
        wd=8-dd;
        for j=1:dif
         z0(wd)=5.D0;
         if (j == dif), break, end
         wd=wd-1;
         if (wd == 0), wd=7; end
        end
       end
       if (Itrad >= 6) 
        for j=1:6
         Y(i,j)=z0(6-j+1)-z0(7);
         if (Ntrad > 0), Y(i,j)=Y(i,j) - Yh(i); end
        end
       elseif ((Itrad == 1) | (Itrad == 2)) 
        sum0=sum(z0(2:6));
        sumf=z0(1)+z0(7);
        if (Ntrad > 0) 
%          sum0=sum0-Z(i,Ntrad1);
%          sumf=sumf+Z(i,Ntrad1);
         sum0=sum0-Yh(i);
         sumf=sumf+Yh(i);
        end
        Y(i,1)=sum0-sumf*double(5)/double(2);
       end
% leap year variable
       if ((Itrad == 7) | (Itrad == 2)) 
        if (z0(8) == 28) 
         Y(i,Itrad)=-.25;
        elseif (z0(8) == 29) 
         Y(i,Itrad)=.75;
        else
         Y(i,Itrad)=0.D0;
        end
       end 
       if (i == N), continue, end
       dd=dd+dif;
       if (dd > 7), dd=dd-7; end
       mm=mm+1;
       if (mm > 12)
        mm=1; yy=yy+1; ly=mod(yy,4);
       end
      end
      if (Mq == 4) 
       N=in;
       for i=1:N
        if (Itrad >= 6) 
         for j=1:6
          ii=(i-1)*3; Y(i,j)=sum(Y(ii+1:ii+3,j));
         end
        elseif ((Itrad == 1) | (Itrad == 2)) 
         ii=(i-1)*3; Y(i,1)=sum(Y(ii+1:ii+3,1));
        end
% leap year variable 
        if ((Itrad == 7) | (Itrad == 2)) 
         ii=(i-1)*3; Y(i,Itrad)=sum(Y(ii+1:ii+3,Itrad));
        end
       end 
       Y(N+1:end,:)=[];
      end
