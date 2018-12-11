% script file to simulate a series that follows a transfer function model.
%The input series is assumed to follow an ARIMA model. The model is
%
%  y_t =            u_t           +               v_t
%      =  (3. + 1.5B)/(1-.5B)x_t  + (1-.4B)(1-.6B^12)/(1-B)(1-B^12)a_t
%
% where x_t follows the model
%
% x_t = (1 - .7B^12)/(1-.5B)(1-B^12)b_t
%
% and std(a_t) = 1. and std(b_t)=.5.
%

clear

freq = 12;

x = arimasimeasy(freq, '[p dr q]', [1, 0, 0], '[ps ds qs]', [0, 1, 1], ...
    'phir', [-.5, 1], 'ths', [-.7, 1], 'N', 150, 'discard', ...
    50, 'seed', 18, 'stda', .5, 'gft', 0);


v = arimasimeasy(freq, '[p dr q]', [0, 1, 1], '[ps ds qs]', [0, 1, 1], ...
    'thr', [-.4, 1], 'ths', [-.6, 1], 'N', 150, 'discard', ...
    50, 'seed', 20, 'stda', 1., 'gft', 0);

%third, filter inputs by Phi(B)^{-1}Gamma(B)
%1) without the input model
thp = [1.5, 3.];
phip = [-.5, 1.];
% u01 = varmafilp(x,phip,thp);
%
% omega(:,:,1)=3.; omega(:,:,2)=1.5;
% delta(:,:,1)=1.; delta(:,:,2)=-.5;
% u02 = varmafilp(x,delta,omega);

%2) with the input model
%set up state space form for input model
phix = [-.5, 1.];
Phix = [-1., 1.];
Thx = [-.7, 1.];
thx = 1.;
Sigma = .5^2;
u1 = varmafilp(x, phip, thp, phix, thx, Phix, Thx, Sigma, freq);

%  phi(:,:,1)=1; phi(:,:,2)=-.5;
%  Phi(:,:,1)=1; Phi(:,:,2)=-1.;
%  Th(:,:,1)=1.; Th(:,:,2)=-.7;
%  th(:,:,1)=1.;
%  u2 = varmafilp(x,delta,omega,phi,th,Phi,Th,Sigma,freq);

u = u1;

y = v + u;

% y(5:9)=NaN(5,1);
% y(30)=NaN;

% out=arimaeasy(x,freq,'pr',1,'gft',1,'sname','myseries');
% out=arimaeasy(v,freq,'pr',1,'gft',1,'sname','myseries');

%Identify and estimate the model
out = tfeasy(y, x, freq, 'pr', 1, 'gft', 1, 'sname', 'mytfseries', 'tfident', 1, ...
    'autmid', 1);