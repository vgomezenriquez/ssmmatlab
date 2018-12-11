function [A,Sigma,Xi]=vincovma(F,K,stda)
%       
%
%        This function computes the elements of the initial state vector
%        for the diffuse Kalman filter in ARIMA models
%
% Copyright (c) 21 July 2003 by Victor Gomez
% Ministerio de Hacienda, Direccion Gral. de Presupuestos, 
% Subdireccion Gral. de Analisis y P.E., 
% Alberto Alcocer 2, 1-P, D-34, 28046, Madrid, SPAIN.
% Phone : +34-915835439
% E-mail: vgomez@sepg.minhac.es
%
% The author assumes no responsibility for errors or damage resulting from the use
% of this code. Usage of this code in applications and/or alterations of it should 
% be referenced. This code may be redistributed if nothing has been added or
% removed and no money is charged. Positive or negative feedback would be appreciated.
%


% compute initial state vector
[Q,T] = schur(F,'real'); 
%eigenvalues are sorted in descending absolute value
[nT,mT]=size(T);
if nT > 1
 [Qt,Tt,ap] = SortSchur(Q,T,0+0i);  
else
 Qt=Q; Tt=T;
end
%extract eigenvalues
if isOctave
 E=sort(eig(Tt));
else
 E=ordeig(Tt); 
end
[n,m]=size(E); nm=max(n,m); 
nu=0;
for i=1:nm
 if abs(E(i)) >= .99d0
  nu=nu+1;
 end
end
% nu
if (nu < nm) & (nu > 0)
 UN=Tt(1:nu,1:nu);  
 US=Tt(nu+1:end,nu+1:end); U12=Tt(1:nu,nu+1:end);
 %el siguiente paso se puede evitar resolviendo la ecuacion 
 %U_NX-XU_S=U_{12} directamente. 
%  X = lyap(UN, -US, -U12); XX=X
 [X,ferror] = mclyapunov(UN,US,U12);
 Q=eye(nm,nm); Q(1:nu,nu+1:end)=X; R=Q*Qt'; 
 PN=Qt(:,1:nu); PS=Qt(:,nu+1:end);
 RN=PN; RS=PS - PN*X;
 A=RN;
 Kt=R*K; KN=Kt(1:nu,:); KS=Kt(nu+1:end,:); GS=KS*stda;
 Sigma = mlyapunov(US,GS*GS',.995); Xi=RS;
elseif (nu == 0)
 G=K*stda;
 A=[]; Sigma = mlyapunov(F,G*G',.995); Xi=eye(nm);
else
 A=F; Xi=[]; Sigma=[];
end
