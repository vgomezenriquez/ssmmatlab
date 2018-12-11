function [str, ferror] = suusm(comp, y, Y, npr)
%************************************************************************
% PURPOSE: this function sets up a univariate structural model given a
% structure containing information about the components. It returns a
% structure containing the model information. The model is
%
%  y_t = Y_t*beta + p_t + s_t + u_t + v_t + e_t,
%
% where Y_t is a vector of regression variables, p_t is the trend, s_t is
% the seassonal, u_t is the cyclical, v_t is the AR, and e_t is the
% irregular component.  This function allows for other values of freq
% besides 1, 4 and 12.
% Note that freq is a field of structures comp and datei. It is also the
% number of seasons. Therefore, it determines the seasonal component.
%---------------------------------------------------
% USAGE: [str,ferror] = suusm(comp,y,Y,npr)
% where structure comp has the following fields:
%          .level  = a 1 x 3 dimensional array such that level(1) is a
%                    code (see below *), level(2) is the standard error of
%                    the level and level(3) = NaN means the standard error
%                    is to be estimated, =0 it is fixed
%          .slope  = a 1 x 3 dimensional array such that slope(1) is a
%                    code (see below *), slope(2) is the standard error of
%                    the slope and slope(3) = NaN means the standard error
%                    is to be estimated, =0 it is fixed
%          .seas   = a 1 x 3 dimensional array such that seas(1) is a
%                    code (see below *), seas(2) is the standard error of
%                    the seasonal and seas(3) = NaN means the standard
%                    error is to be estimated, =0 it is fixed
%          .cycle  = a 1 x 3 dimensional array such that cycle(1) is a code
%                    (see below *), cycle(2) is the standard error of the
%                    cycle and cycle(3) = NaN means the standard error is
%                    to be estimated, =0 it is fixed
%          .cyclep = (only if field .cycle is present) a 2 x 2 array
%                    containing the first row the two cycle parameters (rho
%                    and freqc) and the second row a NaN or zero for each
%                    cycle parameter. Each NaN means that the correspondig
%                    cycle parameter is to be estimated and each zero means
%                    that it is fixed
%          .cycleb = (only if field .cycle is present) a 1 x 2 array
%                    such that cycleb(1) and cycleb(2) contain the end
%                    points of the frequency interval in which the cycle is
%                    supposed to be defined.
%          .ar     = a 1 x 3 dimensional array such that ar(1) is a
%                    code (see below *), ar(2) is the standard error of
%                    the ar component and ar(3) = NaN means the standard
%                    error is to be estimated, =0 it is fixed
%          .arp    = (only if field .ar is present) a 2 x k array, where
%                    k is the order of the autoregressive, containing the
%                    first row the autoregressive parameters and the
%                    second row a NaN or zero for each autoregressive
%                    parameter. Each NaN means that the correspondig
%                    autoregressive parameter is to be estimated and each
%                    zero means that it is fixed
%          .irreg  = a 1 x 3 dimensional array such that irreg(1) is a
%                    code (see below *), irreg(2) is the standard error of
%                    the irregular and irreg(3) = NaN means the standard
%                    error is to be estimated, =0 it is fixed
%          .conout = 'level' if the standard error of the level is
%                    concentrated out
%                    'slope' if the standard error of the slope is
%                    concentrated out
%                    'seas' if the standard error of the seasonal
%                    is concentrated out
%                    'cycle' if the standard error of the cycle is
%                    concentrated out
%                    'ar' if the standard error of the ar component is
%                    concentrated out
%                    'irreg' if the standard error of the irregular is
%                    concentrated out
%                    If .conout is not input, the program will determine
%                    the biggest variance.
%          .freq   = number of observations per year
%          .datei  = calendar structure (output of function cal)
%         .sqrtfil = 0, use ordinary two-stage Kalman filter for estimation
%                    1, use square root version (specially for long series)
%          y       = data vector
%          Y       = matrix for regression variables. It contains the stack
%                    of the Y_t matrices.
%          npr     = number of forecasts
%-------------------------------------------------------------------------
% * codes for the components:
%   level = -1  constant
%            1  stochastic
%            2  Butterworth tangent
%   slope = -1  constant
%            1  stochastic
%   seas  = -1  fixed dummy seasonality
%            1  stochastic dummy seasonality
%            2  trigonometric seasonality
%            4  Butterworth tangent
%   cycle =  1  structural model cycle
%            2  Butterworth sine cycle
%   irreg =  1  stochastic
%   ar    =  k  autoregressive component of order k
%
%---------------------------------------------------
% RETURNS: str = a structure containing the following fields
%
%       matrices according to the model
%
%          y_t = X_t*beta + Z_t*alpha_t + G_t*epsilon_t
%          alpha_{t+1}= W_t*beta + T_t*alpha_t + H_t*epsilon_t
%
%          where epsilon_t is (0,sigma^2I),
%
%          with initial state
%
%          alpha_1= c + W_0*beta + a_1 + A_1*delta
%
%          where c is (0,Omega) and delta is (0,kI) (diffuse)
%      More specifically:
%            .X : an (n x nbeta) matrix containing the X_t matrices;
%                 a  (1 x nbeta) if it is time invariant;
%                 it can be []
%            .Z : an (n x nalpha) matrix containing the Z_t matrices;
%                 a  (1 x nalpha) matrix if it is time invariant
%            .G : an (n x nepsilon) matrix containing the G_t matrices;
%                 a  (1 x nepsilon) matrix if it is time invariant
%            .W : an (n*nalpha x nbeta) matrix containing the W_t matrices;
%                 an (nalpha x nbeta) matrix if it is time invariant;
%                 it can be []
%            .T : an (n*nalpha x nalpha) matrix containing the T_t matrices;
%                 an (nalpha x nalpha) matrix if it time invariant
%            .H : an (n*nalpha x nepsilon) matrix containing the H_t matrices;
%                 an (nalpha x nepsilon) if it is time invariant
%     .ins : an nalpha x (cc+cw0+ca1+cca1) matrix containing the initial
%            state information, according to array i below
%       .i : a  1 x 4 array containing 4 integers, i=[cc cw0 ca1 cca1],
%            where
%            cc   = nalpha if c is not missing (0 if c missing)
%            cw0  = number of columns in W_0 (0 if W_0 missing)
%            ca1  = 1 if a_1 is not missing (0 if a_1 missing)
%            cca1 = number of columns in A_1 (0 if A_1 missing)
%   .trend : trend code
%   .slope : slope code
%    .seas : seasonal code
%   .cycle : cycle code
%     .xl1 : lower bound of the frequency interval in which the
%            cycle is supposed to be defined
%     .xl2 : upper bound of the frequency interval in which the
%            cycle is supposed to be defined
%     .arp : AR code
%   .irreg : irregular code
%    .conc : index for the parameter that is concentrated out
%           (see description below *)
%       .x : vector with all parameters (see description below **)
%      .xv : vector with parameters to be estimated
%      .xf : vector with fixed parameters
%    .pvar : array with variable parameter indices
%    .pfix : array with fixed parameter indices
%   .stord : array containing parameter indices (see description below ***)
%    .freq : frequency of the data
%   .datei : calendar structure
%    .comp : structure comp (input to suusm)
%
%---------------------------------------------------
% * One of the standard deviations is concentrated out and, therefore, is not
%   estimated. The field conout contains information about this standard
%   deviation. The user can select this standard deviation or the program can
%   do it automatically instead. The biggest variance should be selected.
%
% ** we put in x the parameters of the model, except the one that is
%    concentrated out, in the order:
%    1 - irregular standard deviation
%    2 - level standard deviation
%    3 - slope standard deviation
%    4 - seasonal standard deviation
%    5 - autoregressive standard deviation
%    6 - cycle standard deviation
%    7,8 - cycle parameters, rho and frequency
%    9,10,.. autoregressive parameters
%
% *** stord is an index such that its i-th element indicates to which
%     component (according to the ordering above) belongs the i-th element of x.
%*************************************************************************
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
%*************************************************************************

ferror = 0;
str = [];

if isfield(comp, 'datei')
    datei = comp.datei;
else
    disp('field datei should be in structure comp in suusmm');
    ferror = 1;
end

if isfield(comp, 'freq')
    freq = comp.freq;
else
    disp('field freq should be in structure comp in suusmm');
    ferror = 1;
end


% x is split into two subvectors, xv that contains the parameters to be
% estimated, and xf that contains the fixed parameters.

x = [];
xv = [];
xf = [];
pfix = [];
pvar = [];
stord = [];
j = 0;

% Check the validity of entries in fields related to components

% Irregular
if isfield(comp, 'irreg')
    
    if ~isequal(size(comp.irreg), [1, 3])
        disp('field irreg should be an (1 x 3) array')
        ferror = 6;
        return
    end
    
    if any(~isnumeric(comp.irreg))
        disp('all entries in field irreg should be numeric variables')
        ferror = 7;
        return
    end
    
    irregt = comp.irreg(1);
    irreg = comp.irreg(2);
    irregn = comp.irreg(3);
    
    if isnan(irreg)
        disp('irreg(2) has to be a number')
        ferror = 7;
        return
    end
    
    if (irregn ~= 0) && (~isnan(irregn))
        disp('irreg(3) has to be zero or NaN')
        ferror = 7;
        return
    end
    
    if (irregt ~= 1) || ((irregt == 1) && (irreg == 0) && (irregn == 0))
        disp('irregular component can only be stochastic')
        ferror = 4;
        return
    end
    
else
    irregt = 0;
end

% Level
if isfield(comp, 'level')
    
    if ~isequal(size(comp.level), [1, 3])
        disp('field level should be an (1 x 3) array')
        ferror = 6;
        return
    end
    
    if any(~isnumeric(comp.level))
        disp('all entries in field level should be numeric variables')
        ferror = 7;
        return
    end
    
    levelt = comp.level(1);
    level = comp.level(2);
    leveln = comp.level(3);
    
    if (levelt ~= -1) && (levelt ~= 1) && (levelt ~= 2)
        disp('level type should be minus one, one or two')
        ferror = 1;
        return
    end
    
    if isnan(level)
        disp('level(2) has to be a number')
        ferror = 7;
        return
    end
    
    if (leveln ~= 0) && (~isnan(leveln))
        disp('level(3) has to be zero or NaN')
        ferror = 7;
        return
    end
    
    if levelt == -1
        disp('if level(1) = -1, the subsequent entries in comp.level will be ignored')
    end
    
else
    levelt = 0;
end


% Slope
if isfield(comp, 'slope')
    
    if ~isequal(size(comp.slope), [1, 3])
        disp('field slope should be an (1 x 3) array')
        ferror = 6;
        return
    end
    
    if any(~isnumeric(comp.slope))
        disp('all entries in field slope should be numeric variables')
        ferror = 7;
        return
    end
    
    slopet = comp.slope(1);
    slope = comp.slope(2);
    slopen = comp.slope(3);
    
    if (slopet ~= -1) && (slopet ~= 1)
        disp('slope type should be minus one or one')
        ferror = 1;
        return
    end
    
    if isnan(slope)
        disp('slope(2) has to be a number')
        ferror = 7;
        return
    end
    
    if (slopen ~= 0) && (~isnan(slopen))
        disp('slope(3) has to be zero or NaN')
        ferror = 7;
        return
    end
    
    if slopet == 1
        if levelt == -1
            disp('slope(1) cannot be 1 if level(1) = -1. Enter comp.level =  ')
            disp('[1 0 0] if you want to add a slope with the code slope(1) =')
            disp('1 to that level specification')
            ferror = 6;
            return
        end
    elseif slopet == -1
        disp('if slope(1) = -1, the subsequent entries in comp.slope will be ignored')
    end
    
else
    slopet = 0;
end

% Seasonal
if isfield(comp, 'seas')
    
    if ~isequal(size(comp.seas), [1, 3])
        disp('field seas should be an (1 x 3) array')
        ferror = 6;
        return
    end
    
    if any(~isnumeric(comp.seas))
        disp('all entries in field seas should be numeric variables')
        ferror = 7;
        return
    end
    
    seast = comp.seas(1);
    seas = comp.seas(2);
    seasn = comp.seas(3);
    
    if (seast ~= -1) && (seast ~= 1) && (seast ~= 2) && (seast ~= 4)
        disp('seasonal type should be minus one, one, two or four')
        ferror = 1;
        return
    end
    
    if (seast == 4) && (freq ~= 1) && (freq ~= 4) && (freq ~= 12)
        disp('Butterwoth seasonality is only possible with frequencies')
        disp('one, four or twelve')
    end
    
    if isnan(seas)
        disp('seas(2) has to be a number')
        ferror = 7;
        return
    end
    
    if (seasn ~= 0) && (~isnan(seasn))
        disp('seas(3) has to be zero or NaN')
        ferror = 7;
        return
    end
    
else
    seast = 0;
end

% AR
if isfield(comp, 'ar')
    
    if ~isequal(size(comp.ar), [1, 3])
        disp('field ar should be an (1 x 3) array')
        ferror = 6;
        return
    end
    
    if any(~isnumeric(comp.ar))
        disp('all entries in field ar should be numeric variables')
        ferror = 7;
        return
    end
    
    art = comp.ar(1);
    ar = comp.ar(2);
    arn = comp.ar(3);
    
    if isnan(art)
        disp('ar(1) has to be a number')
        ferror = 7;
        return
    elseif (art <= 0)
        disp('AR type should be positive')
        ferror = 1;
        return
    end
    
    if isnan(ar)
        disp('ar(2) has to be a number')
        ferror = 7;
        return
    end
    
    if (arn ~= 0) && (~isnan(arn))
        disp('ar(3) has to be zero or NaN')
        ferror = 7;
        return
    end
    
    if ~isfield(comp, 'arp')
        disp('field arp should be present to specify AR model')
        ferror = 5;
        return
    elseif isfield(comp, 'arp')
        
        if ~isequal(size(comp.arp), [2, art])
            disp('field arp should be an (2 x k) array, where k is given by ar(1)')
            ferror = 6;
            return
        end
        
        if any(~isnumeric(comp.arp))
            disp('all entries in field arp should be numeric variables')
            ferror = 7;
            return
        end
        
        arp = comp.arp(1, :);
        arpn = comp.arp(2, :);
        
        if any(isnan(arp))
            disp('all entries in the first row of field arp should be numbers')
            ferror = 7;
            return
        end
        
        if any(arpn)
            disp('all entries in the second row of field arp should be zero or NaN')
            ferror = 7;
            return
        end
    end
    
elseif ~isfield(comp, 'ar')
    if isfield(comp, 'arp')
        disp('field ar should be present to specify AR model')
        ferror = 5;
    end
    art = 0;
end


% Cycle
if isfield(comp, 'cycle')
    
    if ~isequal(size(comp.cycle), [1, 3])
        disp('field cycle should be an (1 x 3) array')
        ferror = 6;
        return
    end
    
    if any(~isnumeric(comp.cycle))
        disp('all entries in field cycle should be numeric variables')
        ferror = 7;
        return
    end
    
    if ~isfield(comp, 'cyclep')
        disp('field cyclep should be present to specify a cycle')
        ferror = 5;
        return
    elseif isfield(comp, 'cyclep')
        if ~isequal(size(comp.cyclep), [2, 2])
            disp('field cyclep should be an (2 x 2) array')
            ferror = 6;
            return
        end
        if any(~isnumeric(comp.cyclep))
            disp('all entries in field cyclep should be numeric variables')
            ferror = 7;
            return
        end
        
        cyclep = comp.cyclep(1, :);
        cyclepn = comp.cyclep(2, :);
        
        if any(isnan(cyclep))
            disp('both entries in the first row of field cycle should be numbers')
            ferror = 7;
            return
        end
        
        if any(cyclepn)
            disp('entries in the second row of field cycle should be zero or NaN')
            ferror = 7;
            return
        end
    end
    
    if ~isfield(comp, 'cycleb')
        disp('field cycleb should be present to specify a cycle')
        ferror = 5;
        return
    elseif isfield(comp, 'cycleb')
        if ~isequal(size(comp.cycleb), [1, 2])
            disp('field cycleb should be an (1 x 2) array')
            ferror = 6;
            return
        end
        
        cycleb = comp.cycleb;
        
        if any(isnan(cycleb))
            disp('all entries in field cycleb should be numbers')
            ferror = 7;
            return
        end
        
        xl1 = cycleb(1);
        xl2 = cycleb(2);
        
        if (xl1 > xl2)
            disp('first parameter should be smaller than second parameter in field cycleb')
            ferror = 1;
            return
        end
        if (xl1 < 0.) || (xl1 > pi) || (xl2 < 0.) || (xl2 > pi)
            disp('parameters in field cycleb should be between zero and pi')
            ferror = 1;
            return
        end
    end
    
    cyclet = comp.cycle(1);
    cycle = comp.cycle(2);
    cyclen = comp.cycle(3);
    
    if (cyclet ~= 1) && (cyclet ~= 2)
        disp('cycle type should be one or two')
        ferror = 1;
        return
    end
    
elseif ~isfield(comp, 'cycle')
    if (isfield(comp, 'cyclep')) || (isfield(comp, 'cycleb'))
        disp('field cycle should be present to specify a cycle')
        ferror = 5;
    end
    cyclet = 0;
end

if (levelt == -1) && (slopet == -1) && (seast == -1) && (cyclet == 0) ...
        && (art == 0)
    disp('there must be some stochastic component present in the model')
    disp('apart from the irregular. Otherwise, it is a regression model')
    ferror = 6;
    return
end


%********************************************************
% If the user has not selected a variance to be concentrated out, the
% program will find it. The variance to be concentrated out is the biggest
% one.

if ~isfield(comp, 'conout')
    maxconc = -1.d10;
    conc = 0;
    if isfield(comp, 'irreg')
        if isnan(comp.irreg(3)) && (comp.irreg(2) > maxconc)
            maxconc = comp.irreg(2);
            conc = 1;
        end
    end
    if isfield(comp, 'cycle')
        if isnan(comp.cycle(3)) && (comp.cycle(2) > maxconc)
            maxconc = comp.cycle(2);
            conc = 6;
        end
    end
    if isfield(comp, 'ar')
        if isnan(comp.ar(3)) && (comp.ar(2) > maxconc)
            maxconc = comp.ar(2);
            conc = 5;
        end
    end
    if isfield(comp, 'level')
        if isnan(comp.level(3)) && (comp.level(2) > maxconc)
            maxconc = comp.level(2);
            conc = 2;
        end
    end
    if isfield(comp, 'slope')
        if isnan(comp.slope(3)) && (comp.slope(2) > maxconc)
            maxconc = comp.slope(2);
            conc = 3;
        end
    end
    if isfield(comp, 'seas')
        if isnan(comp.seas(3)) && (comp.seas(2) > maxconc)
            %    maxconc=comp.seas(2);
            conc = 4;
        end
    end
    if (conc == 1)
        comp.conout = 'irreg';
    elseif (conc == 2)
        comp.conout = 'level';
    elseif (conc == 3)
        comp.conout = 'slope';
    elseif (conc == 4)
        comp.conout = 'seas';
    elseif (conc == 5)
        comp.conout = 'ar';
    else
        comp.conout = 'cycle';
    end
    conout = comp.conout;
else
    % The user has selected a variance to be concentrated out
    if ~strcmp(comp.conout, {'irreg', 'level', 'slope', 'seas', 'ar', 'cycle'})
        disp('conout should be one of the following strings:')
        disp('irreg, level, slope, seas, ar, cycle')
        ferror = 1;
        return
    end
    conout = comp.conout;
    conc = [];
end

% conout,conc,pause

%************************************************************************

% Put initial parameter values entered by the user into xv, xf and xvf
% Create pvar and pfix

if isfield(comp, 'irreg')
    
    if ~strcmp(conout, 'irreg')
        x = [x, irreg];
        stord = [stord, 1];
        j = j + 1;
        if isnan(irregn)
            xv = [xv, irreg];
            pvar = [pvar, j];
        else
            xf = [xf, irreg];
            pfix = [pfix, j];
        end
    else
        if ~isnan(irregn)
            disp('if irregular s.d. is concentrated out, it cannot be fixed')
            ferror = 3;
            return
        end
        conc = 1;
    end
end


if isfield(comp, 'level')
    
    if levelt > 0
        
        if ~strcmp(conout, 'level')
            x = [x, level];
            stord = [stord, 2];
            j = j + 1;
            if isnan(leveln)
                xv = [xv, level];
                pvar = [pvar, j];
            else
                xf = [xf, level];
                pfix = [pfix, j];
            end
        else
            if ~isnan(leveln)
                disp('if level s.d. is concentrated out, it cannot be fixed')
                ferror = 3;
                return
            end
            conc = 2;
        end
        
    else
        if strcmp(conout, 'level')
            disp('if level s.d. is concentrated out, level(1) cannot be -1')
        end
    end
end


if isfield(comp, 'slope')
    
    if slopet == 1
        if ~strcmp(conout, 'slope')
            x = [x, slope];
            stord = [stord, 3];
            j = j + 1;
            if isnan(slopen)
                xv = [xv, slope];
                pvar = [pvar, j];
            else
                xf = [xf, slope];
                pfix = [pfix, j];
            end
        else
            if ~isnan(slopen)
                disp('if slope s.d. is concentrated out, it cannot be fixed')
                ferror = 3;
                return
            end
            conc = 3;
        end
        
    else
        if (slopet == -1) && strcmp(conout, 'slope')
            disp('if slope s.d. is concentrated out, slope cannot be deterministic')
        end
    end
end


if isfield(comp, 'seas')
    if seast > 0
        if ~strcmp(conout, 'seas')
            x = [x, seas];
            stord = [stord, 4];
            j = j + 1;
            if isnan(seasn)
                xv = [xv, seas];
                pvar = [pvar, j];
            else
                xf = [xf, seas];
                pfix = [pfix, j];
            end
        else
            if ~isnan(seasn)
                disp('if seasonal s.d. is concentrated out, it cannot be fixed')
                ferror = 3;
                return
            end
            conc = 4;
        end
    end
end


if isfield(comp, 'ar')
    
    if ~strcmp(conout, 'ar')
        x = [x, ar];
        stord = [stord, 5];
        j = j + 1;
        if isnan(arn)
            xv = [xv, ar];
            pvar = [pvar, j];
        else
            xf = [xf, ar];
            pfix = [pfix, j];
        end
    else
        if ~isnan(arn)
            disp('if autoregressive s.d. is concentrated out, it cannot be fixed')
            ferror = 3;
            return
        end
        conc = 5;
    end
    
    x = [x, arp];
    narp = length(arp);
    
    %the next line changed 11-03-2012
    %  jj=5+ncycle+ncyclep;
    jj = 8;
    %end of modification
    
    for i = 1:narp
        jj = jj + 1;
        stord = [stord, jj];
        j = j + 1;
        if isnan(arpn(i))
            xv = [xv, arp(i)];
            pvar = [pvar, j];
        else
            xf = [xf, arp(i)];
            pfix = [pfix, j];
        end
    end
end


if isfield(comp, 'cycle')
    
    if ~strcmp(conout, 'cycle')
        x = [x, cycle];
        stord = [stord, 6];
        j = j + 1;
        if isnan(cyclen)
            xv = [xv, cycle];
            pvar = [pvar, j];
        else
            xf = [xf, cycle];
            pfix = [pfix, j];
        end
    else
        if ~isnan(cyclen)
            disp('if the s.d. of the cycle is concentrated out, it cannot be fixed')
            ferror = 3;
            return
        end
        conc = 6;
    end
    
    x = [x, cyclep];
    ncyclep = length(cyclep);
    
    cyclep1 = cyclep(1);
    cyclep2 = cyclep(2);
    if (cyclep1 < 0.) || (cyclep1 > 1.)
        disp('first parameter in field cyclep should be between zero and one')
        ferror = 1;
        return
    end
    
    if (cyclep2 < 0.) || (cyclep2 > pi)
        disp('second parameter in field cyclep should be between zero and pi')
        ferror = 1;
        return
    end
    
    jj = 6;
    
    for i = 1:ncyclep
        jj = jj + 1;
        stord = [stord, jj];
        j = j + 1;
        if isnan(cyclepn(i))
            xv = [xv, cyclep(i)];
            pvar = [pvar, j];
        else
            xf = [xf, cyclep(i)];
            pfix = [pfix, j];
        end
    end
end

% x,xv,xf,pause

%******************************************************************
% all standard deviations should be such that the biggest one is
% equal to one and the rest are equal to or less than 0.1.
stordc = [stord, conc];
nr = length(stordc);
maxabsst = -1.d10;
cc = {'comp.irreg', 'comp.level', 'comp.slope', ...
    'comp.seas', 'comp.ar', 'comp.cycle'};

for i = 1:nr
    if (stordc(i) <= 6)
        stdc = eval(cc{stordc(i)});
        absst = abs(stdc(2));
        if (absst > .1) && (absst > maxabsst)
            maxabsst = absst;
        end
    end
end

if (maxabsst > .1) %change of scale
    for i = 1:nr
        ivarx = 0;
        if stordc(i) == 1
            comp.irreg(2) = comp.irreg(2) * (.1 / maxabsst);
            if conc ~= 1
                ivarx = 1;
            end
        elseif stordc(i) == 2
            comp.level(2) = comp.level(2) * (.1 / maxabsst);
            if conc ~= 2
                ivarx = 1;
            end
        elseif stordc(i) == 3
            comp.slope(2) = comp.slope(2) * (.1 / maxabsst);
            if conc ~= 3
                ivarx = 1;
            end
        elseif stordc(i) == 4
            comp.seas(2) = comp.seas(2) * (.1 / maxabsst);
            if conc ~= 4
                ivarx = 1;
            end
        elseif stordc(i) == 5
            comp.ar(2) = comp.ar(2) * (.1 / maxabsst);
            if conc ~= 5
                ivarx = 1;
            end
        elseif stordc(i) == 6
            comp.cycle(2) = comp.cycle(2) * (.1 / maxabsst);
            if conc ~= 6
                ivarx = 1;
            end
        end
        if ivarx == 1
            x(i) = x(i) * (.1 / maxabsst);
        end
    end
    xv = x(pvar);
    xf = x(pfix);
end
% comp,x,xv,xf,pause

%**********************************************************************
modescr.trend = levelt;
modescr.slope = slopet;
modescr.seas = seast;
modescr.cycle = cyclet;
if cyclet > 0
    modescr.xl1 = xl1;
    modescr.xl2 = xl2;
end
modescr.ar = art;
modescr.irreg = irregt;
modescr.stord = stord;
modescr.conc = conc;

str = modelstruc(xv, y, Y, pfix, pvar, xf, modescr, datei, npr);
str.x = x;
str.xv = xv;
str.xf = xf;
str.pfix = pfix;
str.pvar = pvar;
str.freq = freq;
str.datei = datei;
str.comp = comp;
if isfield(comp, 'sqrtfil')
    str.sqrtfil = comp.sqrtfil;
end