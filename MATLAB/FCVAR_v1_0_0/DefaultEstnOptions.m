function [ FCVARoptions ] = DefaultEstnOptions(db0, k, r, p)
% 
% function [ FCVARoptions ] = DefaultEstnOptions(db0, k, r, p)
% Lee Morin & Morten Nielsen
% August 21, 2011
%
% 
% DefaultEstnOptions(db0, k, r, p) returns a set of preliminary setings
%       corresponding to the default set of options specified below.
% 
% These options are later processed by AdjustEstnOptions(db0, k, r, p, initialFCVARoptions),
%       which completes the specification of the options array.
% 
% 
% input = vector or scalar db0 of fractional differencing parameters.
%       The size of db0 indicates whether or not restriction d = b is imposed
%           based on the number of free parameters.
%       scalars k, r and p denote the number of lags, cointegrating rank
%           and the number of variables, respectively.
% 
% 
% output = cell array FCVARoptions of options for use in estimation routines.
% 
% FCVARoptions = { dbEstnOptions, EstnOptions }, an array of 2 arrays.
% 
% dbEstnOptions is an array of options for the estimation of the fractional differencing parameters d and b.
% dbEstnOptions = { constrained, C_db, c_db, R_db, r_db, dbFminOptions };
%   constrained is an indicator for the restriction that dbMax >= d >= b > dbMin.
%   C_db and c_db are the matrices that specify the above restriction.
%   dbMax is an upper bound for the estimation of d and b.
%   dbMin is a lower bound on d and b to convert a restriction of the form 
%       d >= b >= dbMin to represent the restriction d >= b > 0 with strict inequality.
%   R_db and r_db are optional matrices for defining restrictions of the form R_db*db = r_db.
%       These are left blank in the default model.
%   dbFminOptions is an array of options for built-in matlab functions.
% 
% EstnOptions is an array of options for estimating the remaining parameters.
% EstnOptions = { N, M, deterministics, R_Gamma, r_Gamma, [], [], [],
%                   nCols, A, H, a, b, h, [], [], [], [],
%                   [], [], [], print2screen, printGammas,
%                   restrictFminOptions }
%       is a cell array that carries optional arguments for estimation.
%   scalar N is the length of the matrix of initial values of x.
%   scalar M is the number of MA lags in the fractional differencing operator.
%       If M is positive, the fractional differencing filter is truncated at lag M.
%   string deterministics is a switch to determine if, eg. a 'restricted constant'
%       term should be added to the model.
% 
%   R_Gamma and r_Gamma are optional matrices for defining restrictions of the
%       form R_Gamma*Gamma = r_Gamma (intended for exclusion restrictions only).
%       These are left blank in the default model.
% 
%   nCols, A, H, a, b, h are optional matrices for use when estimating with
%       Johansens's switching algorithm.
%       These are left blank in the default model.
% 
%   print2screen is a scalar indicator to print results (from RankTests and FCVARestn 
%      to screen and to generate graphs of characteristic roots.
%   printGammas is a scalar indicator to print coefficients of lagged differences as well.
% 
%   restrictFminOptions is an array of options for optimization over alpha and beta in the 
%       switching algorithm implemented in GetParamsSwitching.
%______________________________________________________



% Specify the folder containing FCVARtools.
% Add a new directory to the search path in Windows. 
% path(path,'c:/PathToFCVARestn/FCVARtools')
% Add a new directory to the search path in Linux. 
% path(path,'/home/PathToFCVARestn/FCVARtools')

path(path,'C:/FCVAR/FCVARtools')


% Set optimization settings for the first stage optimization (d and b).
% A second stage is estimated later in case [ alpha, beta, or rho ] are restriced.

% Specify options for matlab optimization functions.
dbFminOptions = optimset('MaxFunEvals', 10, 'TolFun', 1e-15, 'TolX', 1e-15, 'Display', 'off');
% dbFminOptions = optimset('MaxFunEvals', 5000, 'TolFun', 1e-6, 'TolX', 1e-6);

% Initial values.
% db0 = [ 1 ];        % One parameter implies restriction d = b.
db0 = [ 1 1 ];      % Two parameters define d and b separately.

% Additional optimization settings.
warning off;
format('long');



% Specify constraints for optimization with fractional differencing parameters.

% Constrain optimization so that d and b satisfy dbMax >= d >= b > 0.
% Note: The tolerance dbMin is introduced to impose an inequality restriction
%   using the equality restriction dbMax >= d >= b >= dbMin.
constrained = 1;
dbMax = 10;
dbMin = 0.01;


% Collect into one array of options for this specification.
dbEstnOptions = { constrained, dbMax, dbMin, [], [], dbFminOptions };

% Later, AdjustEstnOptions function generates C_db and c_db from db0 and constrained
%   which constrains d and b so that C_db*dbHat <= c_db.
% Takes in (Initial) dbEstnOptions, starting with FCVARoptions created here
%   and outputs the same with blanks filled in.
% AdjustEstnOptions returns:
% dbEstnOptions = { constrained, C_db, c_db, R_db, r_db, dbFminOptions };



% Specify additional options for estimation.

% Number of observations reserved for presample period.
N = 50;
% Truncation of MA lags in fractional differencing operator.
M = 0;     % M = 0 is ignored.

% Specify deterministics in the fractional AR model.
deterministics = '         none      ';
% deterministics =   'restricted constant';


% Default model places no restrictions on Gamma.
R_Gamma = [];
r_Gamma = [];


% Default model also places no restrictions on either alpha or beta.
nCols = [];
A = [];
H = [];
a = [];
b = [];
h = [];




% Specify an indicator to print results to screen and generate graphs.
print2screen = 1;
% Set a separate indicator to print lagged coefficients as well.
printGammas = 1;




% Specify options for matlab optimization functions used in restricted estimation.
% restrictFminOptions = optimset('MaxFunEvals', 10000, 'TolFun', 1e-15, 'TolX', 1e-15, 'Display', 'off');
restrictFminOptions = optimset('MaxFunEvals', 1000, 'TolFun', 1e-4, 'TolX', 1e-4, 'Display', 'iter');
% Note: These specifications apply to the optimization performed in GetParamsSwitching().





% Collect these options in a single array.
% After AdjustEstnOptions:
% EstnOptions = { N, M, deterministics, R_Gamma, r_Gamma, [], [], [],
%                   nCols, A, H, a, b, h, Abar, Aperp, replaceCols,
%                   keepCols, [], [], [], print2screen, printGammas, restrictFminOptions },
%     the rest of which is filled in by the function AdjustEstnOptions.
EstnOptions = { N, M, deterministics, R_Gamma, r_Gamma, [], [], [],...
                    nCols, A, H, a, b, h, [], [], [], [],...
                    [], [], [], print2screen, printGammas, restrictFminOptions };



% The options are carried through all functions as needed.
% dbEstnOptions carries the options for the first stage estimation of fractional parameters.
% EstnOptions carries the options for all other functions and restrictions.
FCVARoptions = { dbEstnOptions, EstnOptions };


% Complete the specification of options using AdjustEstnOptions.
[ FCVARoptions ] = AdjustEstnOptions(db0, k, r, p, FCVARoptions);




% end