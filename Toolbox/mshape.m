function SigBar = mshape(Tm)
%This function computes the shape of a matrix

[f, c] = size(Tm);
SigBar = zeros(1, c);
for j = 1:c
    SigBar(j) = find(abs(Tm(:, j)) > 1.e-9, 1, 'first');
end


% function SigBar = mshape(Tm,tol)
% %This function computes the shape of a matrix
%
% if nargin < 2
%  [m,n]=size(Tm);
%  normaa=norm(Tm); tol=double(m*n)*normaa*eps(normaa);
% end
%
%
% [f,c] = size(Tm);
% SigBar=zeros(1,c);
% for j=1:c
%     SigBar(j) = find(abs(Tm(:,j))>tol,1,'first');
% end