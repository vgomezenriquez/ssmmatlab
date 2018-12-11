function ierror = chkstainv(Fs)
% new function: 21-1-2011
% This function checks wether the matrix Fs has eigenvalues with modulus
% greater than or equal to one.
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
[np, mp] = size(Fs);
ierror = 0;
if (np ~= mp)
    disp('wrong dimensions of Fs in chkstainv');
    ierror = 1;
    return
end

eigmod = 1.d0;
[U, Fs] = schur(Fs, 'complex');
nonstat = size(find(abs(diag(Fs)) >= eigmod), 1);
if nonstat > 0
    ierror = 2;
end
% return

%old function:
%function ierror=chkstainv(Phi)
% This function checks whether the matrix polynomial Phi(z) = Phi_0+
% Phi_1z + .... + Phi_pzP^p
% is stable or not.
%
% [np,mp,kp]=size(Phi);
% ierror=0;
% if ( np ~= mp )
%     disp('wrong dimensions of Phi in chkstainv');
%     ierror=1;    return
% end
% kp=2*kp;
% s=np; th=zeros(s,s,1); th(:,:,1)=eye(s);
% Sigma=eye(s);
% [c,ierr]=macgf(Phi,th,Sigma,kp+1);
% [r,ierr]=macrf(c);
% if ierr == 2
%  ierror=2; return
% else
%  for i=1:s
%   for j=1:s
%    if (abs(r(i,j,1)) > 1) & (i ~= j)
%     ierror=3; return
%    end
%   end
%  end
%  for k=2:kp+1
%   for i=1:s
%    for j=1:s
%     if abs(r(i,j,k)) > 1
%      ierror=3; return
%     end
%    end
%   end
%  end
% end