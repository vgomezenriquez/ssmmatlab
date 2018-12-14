function [z,sz] = varmafilp(u,delta,omega,phi,th,Phi,Th,Sigma,freq) 
% This function filters the series u_t  using the filter H(z) =
% omega(z)^(-1)*delta(z), where omega(z)=omega_0 + omega_1*z + omega_2*z^2  
% + ....+omega_q*z^q and delta(z) = I + delta_1*z + ... + delta_p*z^p. The
% series u can follow a VARMApqPQ model or not. If it follows no model, it
% is filtered using zeros as starting values. The series followed by u is
% multiplicative seasonal of the form
%
%   phi(B)Phi(B)u_t = th(B)Th(B)a_t
%
% If there are unit roots, they are assumed to be in phi, Phi or delta
%
% Input arguments:
%     u: the input series
% omega: the filter denominator, as a MATLAB polynomial in the scalar case,
%        and as a matrix polynomial in the vector case
% delta: the filter numerator, as a MATLAB polynomial in the scalar case,
%        and as a matrix polynomial in the vector case
%   phi: the regular AR part, as a MATLAB polynomial in the scalar case,
%        and as a matrix polynomial in the vector case
%   Phi: the seasonal AR part, as a MATLAB polynomial in the scalar case,
%        and as a matrix polynomial in the vector case
%    th: the regular MA part, as a MATLAB polynomial in the scalar case,
%        and as a matrix polynomial in the vector case
%    Th: the seasonal MA part, as a MATLAB polynomial in the scalar case,
%        and as a matrix polynomial in the vector case
% Sigma: the covariance matrix of the input innovations
%  freq: the number of seasons
%
% Output arguments:
% z  : the filtered series
% sz : the mse of z
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

z=[]; sz=[];

if nargin > 3 && nargin < 9
 error('Number of arguments in varmafilp should be 3 or 9')
elseif nargin == 3
 mo=max(length(size(omega)),length(size(delta)));
 if mo == 2
  mdelta=size(delta,2); momega=size(omega,2);
  if mdelta < momega
   delta=[zeros(1,momega-mdelta) delta];
  elseif momega < mdelta
   omega=[zeros(1,mdelta-momega) omega];
  end
  inc= 0;
  [Fs,Bs,Hs]=akaikessm1(delta,omega);
  Ds=omega(end); krog=size(Fs,1);
  % Fs,Bs,Hs,Ds,pause
  z = varmafil(u,Fs,Hs,Bs,Ds,krog,inc);
 else
  m=size(omega,2);
  s=size(omega,1);
  [phige,thetage,krog] = pecheform(delta,omega);
  strg = matechelon(krog,s,m);
  strg = vecparwr(strg); 
  strg.phis=phige; 
  strg.thetas=phige; 
  strg.gammas=thetage;
  strg = armaxe2sse(strg); 
  %Fs, Hs, Bs and Ds are the filter matrices for the inputs
  Fs=strg.Fs; Hs=strg.Hs; Bs=strg.Bs; Ds=strg.Ds;
  inc=0;
  z = varmafil(u,Fs,Hs,Bs,Ds,krog,inc);
 end
elseif nargin == 9
 mo=max(length(size(omega)),length(size(delta))); 
 mi=max([length(size(phi)) length(size(th))...
         length(size(Phi)) length(size(Th))]); 
 if mo ~= mi
  error('filter model and input model should be of the same type')
 end
 if mo == 2
  %first input filter
  mdelta=size(delta,2); momega=size(omega,2);
  if mdelta < momega
   delta=[zeros(1,momega-mdelta) delta];
  elseif momega < mdelta
   omega=[zeros(1,mdelta-momega) omega];
  end
  [Ff,Gf,Hf]=akaikessm1(delta,omega);
  Jf=omega(end); 
  %then, input model
  mphi=size(phi,2);
  mth=size(th,2);
  mPhi=size(Phi,2);
  mTh=size(Th,2);
  ndelta=0;
  if mphi > 1
   arphi=abs(roots(fliplr(phi)));
   for i=1:mphi-1
    if arphi(i) >= .99d0
     ndelta=ndelta+1;
    end
   end
  end
  if mth > 1
   arth=abs(roots(fliplr(th)));
   for i=1:mth-1
    if arth(i) >= .99d0
     ndelta=ndelta+1;
    end
   end
  end
  if mPhi > 1
   arPhi=abs(roots(fliplr(Phi)));
   for i=1:mPhi-1
    if arPhi(i) >= .99d0
     ndelta=ndelta+freq;
    end
   end
   Phic=zeros(1,freq*(mPhi-1)+1);
   Phic(end)=1;
   for i=1:mPhi-1
    Phic((i-1)*freq+1)=Phi(i);
   end
  else
   Phic=1.;
  end
  if mTh > 1
   arTh=abs(roots(fliplr(Th)));
   for i=1:mTh-1
    if arTh(i) >= .99d0
     ndelta=ndelta+freq;
    end
   end
   Thc=zeros(1,freq*(mTh-1)+1);
   Thc(end)=1.;
   for i=1:mTh-1
    Thc((i-1)*freq+1)=Th(i);
   end
  else
   Thc=1.;
  end
  phix=conv(phi,Phic); 
  thetax=conv(th,Thc);
  stda=sqrt(Sigma);
  mphix=size(phix,2); mthetax=size(thetax,2);
  if mphix < mthetax
   phix=[zeros(1,mthetax-mphix) phix];
  elseif mthetax < mphix
   thetax=[zeros(1,mphix-mthetax) thetax];
  end
  [Fx,Gx,Hx]=akaikessm1(phix,thetax);
  Jx=1.;  
  %set up state space form for filtered series
  [Fsp,Gsp,Hsp,Jsp]=cascadessm1(Fx,Gx,Hx,Ff,Gf,Hf,Jf); 
  Gsp=Gsp*stda; Jsp=Jsp*stda;  
  X=[]; W=[];
  T=Fsp; H=Gsp; 
  nx=size(Hx,1); mf=size(Hf,2);
  Z=[Hsp; [zeros(nx,mf)  Hx]]; G=[Jsp; Jx*stda]; 
  %number of unit roots equal to ndelta
  [ins,ii,ferror]=incossm(T,H,ndelta);
  mz=size(Hsp,1);
  uz=[NaN(size(u,1),mz) u];
  %  disp('******** Computation of filtered series   ********');
    %
    % Computation with function smoothgen.m
    % 
    % Function smoothgen smooths a general vector:
    % Y_t = U_t*beta + C_t*alpha_t + D_t*epsilon_t
    % In this case, it is desired to smooth:
    % Y_t = Z*alpha_t + G*epsilon_t
    % mucd = dimension of Y_t
    % Hence, U_t = X, C_t=Z and D_t=G
    U = [];
    mucd=mz;   
    C = Z(1:mz,:);
    D = G(1:mz,:);
    [zint,szint] = smoothgen(uz,X,Z,G,W,T,H,ins,ii,mucd,U,C,D);
    z=zint;
    sz=szint;
 else
  %first input filter
  [Hf,Ff,Gf,Jf,ierrorf] = qarmax2ss2(delta,omega);
  %then, input model
  L=eye(size(Sigma)); str.freq=freq;
  [~,~,Hx,Fx,Gx,Jx,ferror] = varmapqPQ2ssm(phi,th,Phi,Th,L,str);
  %set up state space form for filtered series
  [Fsp,Gsp,Hsp,Jsp]=cascadessm1(Fx,Gx,Hx,Ff,Gf,Hf,Jf);
  stda=chol(Sigma,'lower');
  Gsp=Gsp*stda; Jsp=Jsp*stda;
  X=[]; W=[];
  T=Fsp; H=Gsp; 
  nx=size(Hx,1); mf=size(Hf,2);
  Z=[Hsp; [zeros(nx,mf)  Hx]]; G=[Jsp; Jx*stda]; 
  %number of unit roots equal to nu
  % compute nu
  [Q,TT] = schur(T,'real'); 
  %eigenvalues are sorted in descending absolute value
  nT=size(TT,1);
  if nT > 1
   [~,Tt,~] = SortSchur(Q,TT,0+0i);  
  else
   Tt=TT;
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
  [ins,ii,ferror]=incossm(T,H,nu); 
  mz=size(Hsp,1);
  uz=[NaN(size(u,1),mz) u];
  %  disp('******** Computation of filtered series   ********');
    %
    % Computation with function smoothgen.m
    % 
    % Function smoothgen smooths a general vector:
    % Y_t = U_t*beta + C_t*alpha_t + D_t*epsilon_t
    % In this case, it is desired to smooth:
    % Y_t = Z*alpha_t + G*epsilon_t
    % mucd = dimension of Y_t
    % Hence, U_t = X, C_t=Z and D_t=G
    U = [];
    mucd=mz;   
    C = Z(1:mz,:);
    D = G(1:mz,:);
    [zint,szint] = smoothgen(uz,X,Z,G,W,T,H,ins,ii,mucd,U,C,D);
    z=zint;
    sz=szint;
 end
end


