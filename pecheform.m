function [phie, thetae, kro, ierror] = pecheform(phi, theta, kro)
%
% This function computes the echelon form corresponding to a transfer
% function Psi(z)=phi^{-1}(z)*theta(z) and, possibly, the Kronecker indices.
% It is assumed that phi(z) is square and that phi(0) is nonsingular.
% Polynomial matrix theta(z) can be nonsquare and, therefore, theta(0) is not
% assumed to be the identity matrix.
%---------------------------------------------------
% USAGE: [phie,thetae,kro,ierror] = pecheform(phi,theta,kro)
% where:    phi   = a k x k polynomial matrix with phi(0) nonsingular
%           theta = a k x m polynomial matrix
%           kro   = a 1 x k vector containing the Kronecker indices
%---------------------------------------------------
% RETURNS:
%           phie  = the AR echelon polynomial matrix
%          thetae = the MA echelon polynomial matrix
%           kro   = a 1 x k vector containing the Kronecker indices
%        ierror =1, dimension mismatch in phi and theta
%               =0, there are no errors on input
%---------------------------------------------------
% If kro is not input, the function uses functions housref and nullref on
% the augmented Sylvester matrices constructed with phi and theta to
% compute the Kronecker indices. If kro is input, a system of linear
% equations based on an appropriate augmented Sylvester matrix is solved.
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

phie = [];
thetae = [];
ierror = 0;
[s1, s2, np] = size(phi);
[s3, s4, nq] = size(theta);
if (s1 ~= s2) | (s1 ~= s3)
    disp('wrong dimensions of phi or theta in pecheform');
    ierror = 1;
    return
end
A = phi(:, :, 1); % we enforce phi(0)=I to ensure D(0)=I in the right MFD
if any(any(A-eye(s1)))
    phi(:, :, 1) = eye(s1);
    for i = 2:np
        phi(:, :, i) = A \ phi(:, :, i);
    end
    for i = 1:nq
        theta(:, :, i) = A \ theta(:, :, i);
    end
end

if nargin == 3
    ikro = 1;
else
    ikro = 0;
end
s = s1;
m = s4;

mnpq = max(np, nq) + 2;
if ikro == 0
    %first obtain right MFD
    % [phir,thetar,ierror] = larmax2rarmax(phi,theta,tol); %state space version
    [phir, thetar, kro, ierror] = pleft2rightcmfd(phi, theta, mnpq); %polynomial version
    
    %then, obtain echelon form and Kronecker indices
    [phie, thetae, kro, ierrorb] = pright2leftcmfd(phir, thetar, mnpq);
else
    %first obtain right MFD
    % [phir,thetar,ierror] = larmax2rarmax(phi,theta,tol); %state space version
    [phir, thetar, kror, ierror] = pleft2rightcmfd(phi, theta, mnpq); %polynomial version
    
    %then, obtain echelon form and Kronecker indices
    [phie, thetae, kro, ierrorb] = pright2leftcmfd(phir, thetar, mnpq, kro);
end
% phie = cleanpmat(phie); thetae = cleanpmat(thetae);
return

% %then set basic D and N matrices
% [p1,p2,p3]=size(phir); [t1,t2,t3]=size(thetar);
% d=max(p3,t3); bphir=zeros(p1,p2,d); bthetar=zeros(t1,t2,d);
% bphir(:,:,1:p3)=phir; bthetar(:,:,1:t3)=thetar;
% % check that indeed D(0)=I. If not, change D and N.
% A=bphir(:,:,1);
% if  any(any(A - eye(p1)))
%  bphir(:,:,1)=eye(p1);
%  for i=2:p3
%   bphir(:,:,i)=bphir(:,:,i)/A; %right division because it is a right MFD
%  end
%  for i=1:t3
%   bthetar(:,:,i)=bthetar(:,:,i)/A;
%  end
% end
%
% D=[]; N=[];
% for i=d:-1:1
%  D=[D bphir(:,:,i)]; N=[N bthetar(:,:,i)];
% end
% % d,D,N
%
% Indn=1:s; spm=s+m; %t1=s, p1=m=p2;
% N0=N; DN=D;                          % elements of augmented Sylvester matrix
%
%
% if ikro == 0                         % Kronecker indices are unknown
%  iser=1; cont=0; ckro=0; kro=zeros(1,s); Ili=1:m; % index for independent rows
%  P=zeros(s); T=zeros(s,m); phie(:,:,1)=P; thetae(:,:,1)=T;
%  while iser
%   cont=cont+1;
%   [n0,m0]=size(N0);  N00=N0; Indn0=Indn;
%   Indn=[]; N0=[];
%   cspm=cont*spm;  cspmms=cspm-s;
%   for i=1:n0                          % search last N block from top to bottom
%    DNi=[DN;N00(i,:)];                 % augment DN matrix with row of last N block
%    [Q,R,Indx,ierror] = housref(DNi'); % test whether last row is l.d.
%    if Indx(end) == 1                  % l.d. row found in N block
%     ckro=ckro+1; kro(Indn0(i))=cont-1;% adjust Kronecker indices
%     Idelx=Indn0(i);                   % create index for echelon row
%     [K,ierror] = nullref(R',Indx);    % solve by back substitution to find coefficients
%     KK=zeros(1,cspm);                 % prepare vector with zeros for echelon row
%     lastc=cspmms+Idelx;               % position of the one in phi(0) in echelon row
%     ix=[Ili lastc];                   % index for the columns in the row that are not zero
%     KK(:,ix)=K;                       % row for echelon form
%     for j=1:cont                      % fill in phi and theta rows
%      tpr=KK(:,(j-1)*spm+1:j*spm);
%      phie(Idelx,:,cont-j+1)=tpr(m+1:end);
%      thetae(Idelx,:,cont-j+1)=-tpr(1:m);
%     end
%    else
%     DN=[DN; N00(i,:)];                 % add row in DN matrix
%     Ili=[Ili cspmms+Indn0(i)];         % adjust index for independent rows
%     Indn=[Indn Indn0(i)];              % store remaining l.i. indices only
%     N0=[N0; N00(i,:)];                 % store l.i. rows of last N block only
%    end
%   end
%   if (cont == max(np,nq)) | (ckro == s)  % stop criteria
%    iser=0;
%   end
%   [nd,md]=size(DN);                    % augment Sylvester matrix with l.i. rows in N block
%   DN=[[DN zeros(nd,m)]; [zeros(m,cont*m) D]];
%   [n0,m0]=size(N0);  [ndn,mdn]=size(DN); % N0 contains the remaining l.i. rows in last N block
%   N0=[zeros(n0,mdn-m0) N0];            % adjust size of N0
%   cspm=cont*spm; Ili=[Ili cspm+1:cspm+m];% adjust index for independent rows
%  end
% else                                  % Kronecker indices are known
%  iser=1; cont=0; ckro=0; Ili=1:m;     % index for independent rows
%  P=zeros(s); T=zeros(s,m); phie(:,:,1)=P; thetae(:,:,1)=T;
%  ldrow=zeros(size(kro));
%  for i=1:s                            % find row numbers for l.d. rows
%   ldrow(i)=kro(i)*spm+m+i;
%  end
% %  ldrow
%  while iser
%   cont=cont+1;
%   [n0,m0]=size(N0);  N00=N0; Indn0=Indn;
%   Indn=[]; N0=[];
%   cspm=cont*spm; cspmms=cspm-s;
%   for i=1:n0                          % search last N block from top to bottom
%    DNi=[DN;N00(i,:)];                 % augment DN matrix with row of last N block
%    ipld=cspmms+Indn0(i);              % index for added row
%    if ldrow(Indn0(i)) == ipld         % added row is l.d.
%     [Q,R,Indx,ierror] = housref(DNi');
%     ckro=ckro+1;                      % adjust Kronecker index count
%     Idelx=Indn0(i);                   % create index for echelon row
%     [K,ierror] = nullref(R',Indx);    % solve by back substitution to find coefficients
%     KK=zeros(1,cspm);                 % prepare vector with zeros for echelon row
%     lastc=cspmms+Idelx;               % position of the one in phi(0) in echelon row
%     ix=[Ili lastc];                   % index for the columns in the row that are not zero
%     KK(:,ix)=K;                       % row for echelon form
%     for j=1:cont                      % fill in phi and theta rows
%      tpr=KK(:,(j-1)*spm+1:j*spm);
%      phie(Idelx,:,cont-j+1)=tpr(m+1:end);
%      thetae(Idelx,:,cont-j+1)=-tpr(1:m);
%     end
%    else
%     DN=[DN; N00(i,:)];                 % add row in DN matrix
%     Ili=[Ili ipld];                    % adjust index for independent rows
%     Indn=[Indn Indn0(i)];              % store remaining l.i. indices only
%     N0=[N0; N00(i,:)];                 % store l.i. rows of last N block only
%    end
%   end
%   if (cont == max(np,nq)) | (ckro == s)  % stop criteria
%    iser=0;
%   end
%   [nd,md]=size(DN);                    % augment Sylvester matrix with l.i. rows in N block
%   DN=[[DN zeros(nd,m)]; [zeros(m,cont*m) D]];
%   [n0,m0]=size(N0);  [ndn,mdn]=size(DN); % N0 contains the remaining l.i. rows in last N block
%   N0=[zeros(n0,mdn-m0) N0];            % adjust size of N0
%   Ili=[Ili cspm+1:cspm+m];             % adjust index for independent rows
%  end
% end
%
