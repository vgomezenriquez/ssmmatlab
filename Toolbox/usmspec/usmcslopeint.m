function ser = usmcslopeint
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
% Series is German consumer real wage series, seasonally adjusted. The
% series is assumed to have a slope intervention in the first quarter of
% 2003 (an impulse).

data = load(fullfile('data', 'PROJECTDATA.dat'));
data(any(isnan(data)'), :) = [];

yor = data(:, 2);
ly = length(yor);
Y = []; %matrix for regression variables
npr = 10; %number of forecasts
freq = 4; % quarterly data
bg_year = 1970;
bg_per = 1;

% Determine the observation number at the time of the intervention
i_year = 2003; % starting year of the intervention
i_per = 1; % starting period of the intervention
% observation number at the time of the intervention
in = ((i_year + i_per) - (bg_year + bg_per)) * freq + 1;

% Incorporate the structural break into the equation for trend slope
% Intervention variable is in this case a pulse variable, i.e. it
% takes the value 1 at the time point of the intervention and 0 otherwise.
%
% The state space model is given by:
% alpha_{t+1} = W_t * beta + T_t * alpha_t + H_t * eps_t
%         Y_t = X_t * beta + Z_t * alpha_t + G_t * eps_t
%
% The univariate structural time series model here is given by:
%     y_t = p_t + u_t + e_t
% p_{t+1} = p_t + b_t + c_t
% b_{t+1} = b_t + d_t
% u_t: trigonometric cycle
%
% The structural change in the trend occurs at tau = 2003.1,
% but the change in slope already occurs at tau-1 = 2002.4:
%    p_{tau-1} =     p_{tau-2} + b_{tau-2} + c_{tau-2}
%    b_{tau-1} =         1 * w + b_{tau-2} + d_{tau-2}
%
%    p_tau =     p_{tau-1} + b_{tau-1} + c_{tau-1}
%    b_tau =                 b_{tau-1} + d_{tau-1}
%
% Given the structure of the system matrices W_t in the state space model,
% it follows that:
%   W_{tau-2} = [0 1 0 0]', where tau is the time point of the intervention
%         W_t = [0 0 0 0]' for t ~= tau-2
% Construct super matrix W consisting of the matrices W_t

W = zeros(ly*4, 1);
inW = in - 2;
W(((inW - 1) * 4 + 1):(((inW - 1) * 4) + 1 + 3)) = [0, 1, 0, 0]';


% Specify components, initial values and fixed parameters
comp.level = [1, 0, 0];
comp.slope = [1, 0, 0];
comp.irreg = [1, .1, NaN];
comp.cycle = [1, .1, NaN];
twopi = 2 * pi;
comp.cyclep = [0.9, twopi / 40.; NaN, NaN];
comp.cycleb = [twopi / 60., twopi / 6.];

ser.yor = yor;
ser.Y = Y;
ser.W = W;
ser.bg_year = bg_year;
ser.bg_per = bg_per;
ser.freq = freq;
ser.comp = comp;
ser.npr = npr;
ser.lam = 1;
ser.olsres = 1;
ser.gft = 1;
