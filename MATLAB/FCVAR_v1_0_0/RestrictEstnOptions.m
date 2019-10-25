function [ restrictFCVARoptions ] = RestrictEstnOptions(db0, k, r, p, defaultFCVARoptions)
% 
% function [ restrictFCVARoptions ] = RestrictEstnOptions(db0, k, r, p, defaultFCVARoptions)
% Lee Morin & Morten Nielsen
% August 21, 2011
%
% 
% RestrictEstnOptions(db0, k, r, p, defaultFCVARoptions) returns a set of preliminary 
%       setings corresponding to the restrictions on parameters specified below.
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
%       array of arrays defaultFCVARoptions contains the preliminary
%           specification of FCVARoptions output from DefaultEstnOptions(db0, k, r, p).
% 
% 
% output = cell array restrictFCVARoptions of options for use in estimation routines.
% 
% restrictFCVARoptions = { dbEstnOptions, EstnOptions }, an array of 2 arrays.
% 
% dbEstnOptions is an array of options for the estimation of the fractional differencing parameters d and b.
% dbEstnOptions = { constrained, C_db, c_db, R_db, r_db, dbFminOptions };
%   constrained is an indicator for the restriction that dbMax >= d >= b > dbMin.
%   C_db and c_db are the matrices that specify the above restriction.
%   dbMax is an upper bound for the estimation of d and b.
%   dbMin is a lower bound on d and b to convert a restriction of the form 
%       d >= b >= dbMin to represent the restriction d >= b > 0 with strict
%       inequality.
%   R_db and r_db are optional matrices for defining restrictions of the form R_db*db = r_db.
%   dbFminOptions is an array of options for built-in matlab functions.
% 
% EstnOptions is an array of options for estimating the remaining parameters.
% EstnOptions = { N, M, deterministics, R_Gamma, r_Gamma, [], [], [],
%                   nCols, A, H, a, b, h, Abar, Aperp, replaceCols, keepCols,
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
% 
%   nCols, A, H, a, b, h are optional matrices for use when estimating with Johansens's switching algorithm.
%       nCols is a vector of dimensions r_i of the partitioned matrices of alpha and beta.
%       where sum(r_i) = r.
%   A is an array of linear restrictions on each alpha such that alpha_i = A_i*psi_i.
%   H is an array of linear restrictions on each beta such that beta_i = H_i*phi_i.
%   a is an array of values of alpha considered known.
%   b is an array of values of beta considered known.
%   h is an array of vectors specifying linear restrictions on rho.
% 
%   print2screen is a scalar indicator to print results (from RankTests and FCVARestn 
%      to screen and to generate graphs of characteristic roots.
%   printGammas is a scalar indicator to print coefficients of lagged differences as well.
% 
%   restrictFminOptions is an array of options for optimization over alpha and beta in the 
%       switching algorithm implemented in GetParamsSwitching.
%______________________________________________________



% Extract array of options from unrestricted model.
dbEstnOptions = defaultFCVARoptions{1};
EstnOptions = defaultFCVARoptions{2};



% Specify constraints for optimization with fractional differencing parameters.
% Specify equality restrictions on d and b: R_db*dbHat = r_db.
R_db = [];
r_db = [];
% Example: d = 0.7 (for case of d > b).
% R_db = [ 1 0 ];
% r_db = 0.7;


% Additional restrictions can be set to impose the restriction.
C_db = dbEstnOptions{2};
c_db = dbEstnOptions{3};

% Modify C_db and c_db as desired here. If C_db is a matrix, these matrices
% are not modified in AdjustEstnOptions.
% C_db = ...
% c_db = ...

dbEstnOptions{2} = C_db;
dbEstnOptions{3} = c_db;


% Collect into one array of options for this specification.
% dbEstnOptions = { constrained, C_db, c_db, R_db, r_db, dbFminOptions };
dbEstnOptions{4} = R_db;
dbEstnOptions{5} = r_db;




% Specify additional options in case Gamma or [ alpha, beta, rho ] are restriced.

% Set restrictions on Gamma (pk x p) matrix if required.

% Example1: Restrict the top row of all Gamma matrices:
% R_Gamma = [ eye(p*k) zeros(p*k, (p-1)*p*k) ];
% Example2: Restrict the first 2 rows of the first Gamma matrix:
% R_Gamma = [ eye(p) zeros(p, p*p*k-p);...
%             zeros(p,p*k) eye(p) zeros(p, (p-1)*p*k-p) ];

% Note that only zero restrictions are permitted, otherwise you should
%   estimate Gamma jointly with the other parameters.
% Note that the dimensions must conform to the data.
% Otherwise, if no restrictions, leave R_Gamma and r_Gamma blank.
R_Gamma = [];

% r_Gamma left blank here and created in next function.
% Does not allow user to estimate nonzero restrictions.




% Specify restrictions designed for switching algorithm.


% Specify as an array of matrices, { nCols, A, H, a, b, h }, each of which have blanks
%   wherever not specified.
% Later program will fill in zeros and identities where appropriate
%   and produce Abar and APerp, etc as needed.

% First set list of column widths for partition of alpha and beta matrices.
% Note that each of { nCols, A, H, a, b } must have the same number of elements.
nCols = [ 1 1 ];
% Elements of nCols correspond to r_i's below
% No restrictions implied here if nCols is left blank.



% Choose restrictions for beta compatible with switching algorithm.
% Hypotheses of the form beta = (H_1*phi_1,...H_rH*phi_rH)
%   where H_i is p x s_i, phi_i is s_i x r_i and sum(r_i) = r.

H1 = [ 1 0; 0 0; 1 -1 ];
H2 = [];
% Could make restrictions dependent on p.
% H1 = [ ones(p,1)            [ ones(1,p-2); -eye(p-2); zeros(1,p-2) ] ];
% H2 = [ [ ones(p-1,1); 0 ]   [ ones(1,p-2); -eye(p-2); zeros(1,p-2) ] ];

H = { H1, H2 };






% Choose restrictions for alpha compatible with switching algorithm.
% Hypotheses of the form alpha = (A_1*psi_1, ... , A_rA*psi_rA)
%   where A_i is p x m_i, psi_i is m_i x r_i and sum(r_i) = r.
% Note that alpha partition must conform to that for beta (same r_i's).

A1 = [];
A2 = [ 1 0; 1 1; 1 2 ];
A = { A1, A2 };


% Choose restrictions for known parameters.
% a and b allow for a row of alpha or a column of beta to be fixed and considered known.
% A function will later organize the a, b, A and H matrices (along with h for rho).
% User must take care not to jointly specify linear and known restrictions
% on the same columns or rows of beta or alpha.

a1 = 0.5*[ 1 1 1 ]'; % a1 is known.
a2 = []; % a2 is estimated unrestrictedly OR restrictions are specified elsewhere.
a = { a1, a2 };

b1 = []; % b1 is estimated unrestrictedly OR restrictions are specified elsewhere.
b2 = [ 0 1 1 ]'; % b2 is known.
b = { b1, b2 };


% Choose restrictions on rho, if applicable.
% Note: Ignored if not estimating restricted constants model.
h = [ 1 0 ];
% Only zero-one restrictions apply here.





% Collect these options in a single array.
% EstnOptions = { N, M, deterministics, R_Gamma, r_Gamma, [], [], [],...
%                     nCols, A, H, a, b, h, [], [], [], [],...
%                     uGridMax, uGridStep, charPolyDegree, print2screen, printGammas, restrictFminOptions };
%     the rest of which is filled in by the function AdjustEstnOptions.
EstnOptions{4} = R_Gamma;
EstnOptions{5} = [];
EstnOptions{9} = nCols;
EstnOptions{10} = A;
EstnOptions{11} = H;
EstnOptions{12} = a;
EstnOptions{13} = b;
EstnOptions{14} = h;


% The options are carried through all functions as needed.
% dbEstnOptions carries the options for the first stage estimation of fractional parameters.
% EstnOptions carries the options for all other functions and restrictions.
restrictFCVARoptions = { dbEstnOptions, EstnOptions };


% Complete the specification of options using AdjustEstnOptions.
[ restrictFCVARoptions ] = AdjustEstnOptions(db0, k, r, p, restrictFCVARoptions);


% end