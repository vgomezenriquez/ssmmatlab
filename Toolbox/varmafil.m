function [z,rx1] = varmafil(u,F,H,B,D,kro,inc) 
% This function filters the series u_t  using the filter H(z) =
% omega(z)^(-1)*delta(z), where omega(z)=omega_0 + omega_1*z + omega_2*z^2 +
% ....+omega_q*z^q and delta(z) = I + delta_1*z + ... + delta_p*z^p. 
%
% Input arguments:
% u: the input series
% F: a matrix
% H: a matrix
% B: a matrix
% D: a matrix
% kro: a vector containing the Kronecker indices of the filter H(z)
% such that
% 
%  x_{t+1}} = F*x_{t} + B*u_{t}
%  z_{t}    = H*x_{t} + D*u_{t}
%
% inc: = 1 initial conditions, x_{1}, for the filter are estimated
%      = 0 initial conditions equal to zero
%
% Output arguments:
% z  : the output series
% rx1: a matrix containing the design matrix to estimate x_1
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

z=[]; rx1=[];

so=sum(kro); %system order (McMillan degree)
F=sparse(F); H=sparse(H); D=sparse(D); B=sparse(B);
[nd,md]=size(D);
[nu,mu]=size(u);
z=zeros(nu,nd);

if inc == 1
 %this is to estimate x_1=beta by regression later. rx1 contains the 
 %design matrix for beta.
 A=eye(so); rx1=mulHA(H,A,kro);
 for i=1:nu-1
  A=mulFA(F,A,kro);         
  rx1=[rx1;mulHA(H,A,kro)];  
 end
end


x=zeros(so,1);
%z is the filtered series with x_1=0;
for i=1:nu
 U=u(i,:)';
 V=mulHA(H,x,kro) + D*U;
 z(i,:)=V';
 x=mulFA(F,x,kro) + B*U;
end

