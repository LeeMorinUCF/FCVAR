function [ df, dfWarning ] = DegreesOfFreedom(db0, k, r, p, FCVARoptions)

% function [ df, dfWarning ] = DegreesOfFreedom(db0, k, r, p, FCVARoptions)
% Lee Morin & Morten Nielsen
% August 21, 2011
% 
% DegreesOfFreedom(db0, k, r, p, FCVARoptions) calculates the degrees of
%       freedom for the fractional cointegrated VAR model.
% 
% input = scalar or vector dbParam provides the fractional differencing parameters
%           d and b (d = b if dbParam is a scalar).
%   scalar k is the order of lag polynomial in delta d * Lbk * x.
%   scalar r is the cointegrating rank.
%   scalar p is the dimension of the system.
% 
% 
% FCVARoptions = { dbEstnOptions, EstnOptions }, an array of 2 arrays.
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
% 
% 
% output = scalar df is the degrees of freedom for the fractional cointegrated VAR model.
%   scalar dfWarning is a binary variable to indicate whether df does not include 
%       the degrees of freedom in Pi implied by the parameters alpha and beta. 
% 
% No function dependencies.
% 
% Note: When no restrictions are specified on alpha or beta, the degrees of freedom 
%   are adjusted to restrict the top rows of beta to equal the identity matrix. 
%   For restricted estimation, it is assumed that the user imposes enough restrictions 
%   such that the parameters are identified. 
% 
% Note: Whenever there are restrictions on both alpha and beta, the warning dfWarning is given 
%   to indicate that the output df does not include the number of degrees of freedom in Pi.
%   This happens when there are restrictions placed on both alpha and beta together, 
%   in which case it is possible that some restrictions are redundant. 
% 
% Need to fix ~isempty switches because some restrictions are identities when empty.
%______________________________________________________



% Determine the restrictions on the cointegrating relations.

% Options for estimating fractional differencing parameters.
dbEstnOptions = FCVARoptions{1};
r_db = dbEstnOptions{5};

% Options for model specification and restrictions.
EstnOptions = FCVARoptions{2};
% EstnOptions = { N, M, deterministics, R_Gamma, r_Gamma, [], [], [], 
%                   nCols, A, H, a, b, h, RestFminOptions };
deterministics = EstnOptions(3);
R_Gamma = EstnOptions{4};
r_Gamma = EstnOptions{5};
nCols = EstnOptions{9};
nRestns = length(nCols);
A = EstnOptions{10};
H = EstnOptions{11};
a = EstnOptions{12};
b = EstnOptions{13};
h = EstnOptions{14};



% For the unrestricted model.
% df = length(db0) + p*p*k + p*r + r*(p - r) + r*strcmp(deterministics, 'restricted constant');
% In case of restrictions, df is calculated for each parameter in turn.

dfdb0 = length(db0);

dfGamma = p*p*k;

dfBeta = p*r;

dfAlpha = p*r;

dfRho = r*strcmp(deterministics, 'restricted constant');



% Adjust df for restrictions or ID condition.

dfdb0 = dfdb0 - length(r_db);

if ~isempty(R_Gamma)
    % Restrictions on Gamma are represented by ones in the selection matrix R_Gamma.
    dfGamma = dfGamma - sum(sum( R_Gamma ));
end


% Unrestricted case.
if ( isempty(H) && isempty(A) && isempty(a) && isempty(b) && isempty(h) )
    % Impose the ID condition on beta in case of no restrictions on Pi.
    dfBeta = dfBeta - r*r;    
end


if (1) % Condition for dfWarning.
    
    dfBeta = 0;
    dfAlpha = 0;
    dfWarning = 1;
    
else

    for i = 1:nRestns

        % Adjust alpha, beta, and gamma here.
        Ai = A(i);
        if ~isempty(Ai)
            dfAlpha = dfAlpha - p*nCols(i) + size(Ai, 2)*nCols(i);
        end

        Hi = H(i);
        if ~isempty(Hi)
            dfBeta = dfBeta - size(Hi, 1)*nCols(i) + size(Hi, 2)*nCols(i);
            % Note: p could be replaced with (p+1) if there are restrictions on
            %   rho. In this case, H is actually Hstar, where betaStar = Hstar*phi.
        end

        ai = a(i);
        if ~isempty(ai)
            dfAlpha = dfAlpha - p*nCols(i);
        end

        bi = b(i);
        if ~isempty(bi)
            dfBeta = dfBeta - p*nCols(i);
        end

        hi = h(i);
        if ~isempty(hi)
            dfRho = dfRho - (hi == 0);
            % if ~isempty(Hi)
            % if ~isempty(bi)
        end

    end
    
    dfWarning = 0;
end


df = dfdb0 + dfGamma + dfAlpha + dfBeta + dfRho;



% end