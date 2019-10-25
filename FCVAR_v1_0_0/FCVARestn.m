function  [ estimates ] = FCVARestn(x, k, r, db0, FCVARoptions)

% function  [ estimates ] = FCVARestn(x, k, r, db0, FCVARoptions)
% Lee Morin
% August 21, 2011
%
% Estimates the fractional cointegration model under desired options 
%   and prints results to screen as required.
% 
% input = vector or matrix x of data.
%       vector or scalar db0 of starting values for d and b for the fractional differencing parameters.
%           The size of db0 indicates whether or not restriction d = b is imposed
%           based on the number of free parameters.
%       scalar k denotes the number of lags.
%       scalar r denotes the cointegrating rank.
% 
% output = cell array estimates = { dbHat, alphaHat, betaHat, rhoHat, PiHat,
%                                   OmegaHat, GammaHat }
% 
% FCVARoptions = { dbEstnOptions, EstnOptions }, an array of 2 arrays of options 
%       for use in estimation routines.
% dbEstnOptions is an array of options for the estimation of the fractional differencing parameters d and b.
% dbEstnOptions = { constrained, C_db, c_db, R_db, r_db, dbFminOptions };
% EstnOptions is an array of options for estimating the remaining parameters.
%   EstnOptions = { N, M, deterministics, R_Gamma, r_Gamma, [], [], [],
%                   nCols, A, H, a, b, h, [], [], [], [],
%                   [], [], [], print2screen, printGammas,
%                   restrictFminOptions }
% The elements of EstnOptions are described in detail in the comments in
%       functions DefaultEstnOptions(db0, k, r, p),
%       and RestrictEstnOptions(db0, k, r, p, defaultFCVARoptions)
%       and AdjustEstnOptions(db0, k, r, p, initialFCVARoptions).
% 
% 
% Calls the functions FCVARlike and GetParams for estimation, 
%   the functions FCVARhess, SEmat2vec, and SEvec2mat for calculation of standard errors, 
%   and the functions DegreesOfFreedom and CharPolyRoots for display of the estimation results.
%______________________________________________________





% Determine (initial) dimensions of system.
T = size(x, 1); % Length of sample (before truncation for initial values).
p = size(x, 2); % Number of variables.


% Options for estimating fractional differencing parameters.
dbEstnOptions = FCVARoptions{1};
C_db = dbEstnOptions{2};
c_db = dbEstnOptions{3};
R_db = dbEstnOptions{4};
r_db = dbEstnOptions{5};
dbFminOptions = dbEstnOptions{end};
% Options for model specification and restrictions.
EstnOptions = FCVARoptions{2};
N = EstnOptions{1};
T = T - N;
M = EstnOptions{2};
deterministics = EstnOptions{3};
% Indicator for printing results to screen.
print2screen = EstnOptions{22};
printGammas = EstnOptions{23};

% Compute number of free parameters.
[ df, dfWarning ] = DegreesOfFreedom(db0, k, r, p, FCVARoptions);
if dfWarning
    fprintf(1, '\n\nWarning: Degrees of freedom do not reflect the parameters in PiHat.');
    fprintf(1, '\n Calculate degrees of freedom for PiHat manually and add to df reported below.\n\n');
    fprintf(1, '-----------------------------------------------------------------------------------------------------\n');
end

% Perform estimation.

[ dbHat, maxLike, exitflag ] ...
        = fmincon(@( param ) -FCVARlike(x, k, r, param, EstnOptions), db0,...
            C_db, c_db, R_db, r_db, [], [], [], dbFminOptions );


maxLike = - maxLike;

% Load cell array of parameter estimates.
% estimates = { dbHat, alphaHat, betaHat, rhoHat, PiHat, OmegaHat, GammaHat }. 
[ estimates ] = GetParams(x, k, r, dbHat, EstnOptions);
dbHat = cell2mat(estimates(1));
alphaHat = cell2mat(estimates(2));
betaHat = cell2mat(estimates(3));
rhoHat = cell2mat(estimates(4));
PiHat = cell2mat(estimates(5));
OmegaHat = cell2mat(estimates(6));
GammaHat = cell2mat(estimates(7));

% Think about modifying to deal with restrictions.

% Calculate the standard errors of db, alpha and Gamma.
% Note: No standard errors when d = b binds.
if ( length(dbHat) == 2 ) && ( ( dbHat(1) - dbHat(2) ) < 10^(-6) )
    fprintf(1, '\n\nWarning: Estimates are on the d = b boundary so no standard errors are provided.');
    fprintf(1,   '\n         Re-estimate the model with d = b imposed.\n\n');
    SE = zeros(size( SEmat2vec( dbHat, alphaHat, GammaHat, k, r, p ) ));
    
elseif ( length(dbHat) == 2 ) && ( dbHat(2) < 10^(-6) )
    fprintf(1, '\n\nWarning: Estimates are near the b = 0 boundary so no standard errors are provided.');
    SE = zeros(size( SEmat2vec( dbHat, alphaHat, GammaHat, k, r, p ) ));
    
else
    [ hessian ] = FCVARhess(x, k, r, dbHat, alphaHat, betaHat, rhoHat, GammaHat, EstnOptions);
    SE = sqrt(diag(-inv(hessian)));
end
[ dbSE, alphaSE, GammaSE ] = SEvec2mat( SE, k, r, p );


% Print the results to screen.
if print2screen

    fprintf(1,'\n\n-----------------------------------------------------------------------------------------------------\n');
    fprintf(1,'\n\n                      Fractionally Cointegrated VAR: Estimation Results                              \n');
    fprintf(1,'\n\n-----------------------------------------------------------------------------------------------------\n\n');
    fprintf(1,'Dimension of System:  %8.0f     \t Number of observations for estimation: \t %8.0f \n', p, T);
    fprintf(1,'Deterministics:   %s \t Initial values:                          \t %8.0f\n', char(deterministics), N );
    fprintf(1,'Number of Lags:       %8.0f     \t Number of observations in sample:      \t %8.0f \n', k, T+N);
    if M > 0
        fprintf(1,'\n\nTruncation of fractional filter after %d lags.\n', M);
    end
    fprintf(1,'\n-----------------------------------------------------------------------------------------------------\n\n');
    fprintf(1,'Cointegrating rank:    %8.0f           \t\t AIC:                    \t     %8.3f \n', r, maxLike - df);
    fprintf(1,'Log-likelihood:            %8.3f           \t\t BIC:                    \t     %8.3f \n', maxLike, maxLike - df*log(T)/2);
    fprintf(1,'log(det(Omega_hat)):       %8.3f           \t\t Number of parameters:   \t %8.0f \n', log(det(OmegaHat)), df);
    fprintf(1,'\n-----------------------------------------------------------------------------------------------------\n');
    fprintf(1,    '    Fractional Parameters:                                                                             \n');
    fprintf(1,    '-----------------------------------------------------------------------------------------------------\n');
    fprintf(1,    '    Coefficient              \t Estimate              \t  Standard Error \n');
    fprintf(1,    '-----------------------------------------------------------------------------------------------------\n');
    fprintf(1,    '         d                   \t %8.3f              \t     %8.3f                \n', dbHat(1), dbSE(1));
    if length(dbHat) > 1
        fprintf(1,'         b                   \t %8.3f              \t     %8.3f                \n', dbHat(2), dbSE(2));
    end
    fprintf(1,  '-----------------------------------------------------------------------------------------------------\n');


    % Note: No standard errors when d = b binds.
    if ( length(dbHat) == 2 ) && ( ( dbHat(1) - dbHat(2) ) < 10^(-6) )
        fprintf(1, '\nWarning: Estimates are on the d = b boundary so no standard errors are provided.');
        fprintf(1, '\n         Re-estimate the model with d = b imposed.\n\n');
        fprintf(1, '-----------------------------------------------------------------------------------------------------\n');
    elseif ( length(dbHat) == 2 ) && ( dbHat(2) < 10^(-6) )
        fprintf(1, '\nWarning: Estimates are near the b = 0 boundary so no standard errors are provided.\n\n');
        fprintf(1, '-----------------------------------------------------------------------------------------------------\n');
    end

    if r > 0
        if strcmp(deterministics, 'restricted constant')
            varList = '(beta and rho):';
        else
            varList = '(beta):        ';
        end
        fprintf(1,'\n    Cointegrating Equations %s                                                          \n', varList);
        fprintf(1,  '-----------------------------------------------------------------------------------------------------\n');
        fprintf(1,    '      Variable      ' );
        for j = 1:r
            fprintf(1,    '  CI equation %d  ', j);
        end
        fprintf(1,'\n');
        fprintf(1,  '-----------------------------------------------------------------------------------------------------\n');
        for i = 1:p
            fprintf(1,    '        Var%d       ',i );
            for j = 1:r
                fprintf(1,'    %8.3f     ', betaHat(i,j) );
            end
            fprintf(1,'\n');
        end
        if strcmp(deterministics, 'restricted constant')
            fprintf(1,    '      Constant     ' );
            for j = 1:r
                fprintf(1,'    %8.3f     ', rhoHat(j) );
            end
            fprintf(1,'\n');
        end
        fprintf(1,  '-----------------------------------------------------------------------------------------------------\n');
        fprintf(1,  'Note: Identifying restriction imposed.                                                               \n');
        fprintf(1,  '-----------------------------------------------------------------------------------------------------\n\n');

        fprintf(1,'    Adjustment Matrix (alpha):                                                                         \n' );
        fprintf(1,  '-----------------------------------------------------------------------------------------------------\n');
        fprintf(1,    '      Variable      ' );
        for j = 1:r
            fprintf(1,    '  CI equation %d  ', j);
        end
        fprintf(1,'\n');
        fprintf(1,  '-----------------------------------------------------------------------------------------------------\n');
        for i = 1:p
            fprintf(1,    '        Var %d      ',i );
            for j = 1:r
                fprintf(1,'    %8.3f     ', alphaHat(i,j) );
            end
            fprintf(1,'\n');
            fprintf(1,    '         SE %d      ',i );
            for j = 1:r
                fprintf(1,'   (%8.3f  )  ', alphaSE(i,j) );
            end
            fprintf(1,'\n');
        end
        fprintf(1,  '-----------------------------------------------------------------------------------------------------\n');
        fprintf(1,  'Note: Standard Errors in parenthesis.                                                                \n');
        fprintf(1,  '-----------------------------------------------------------------------------------------------------\n\n');

        fprintf(1,'    Long Run Matrix (Pi):                                                                       \n' );
        fprintf(1,  '-----------------------------------------------------------------------------------------------------\n');
        fprintf(1,    '      Variable  ' );
        for j = 1:p
            fprintf(1,    '       Var %d   ', j);
        end
        fprintf(1,'\n');
        fprintf(1,  '-----------------------------------------------------------------------------------------------------\n');
        for i = 1:p
            fprintf(1,    '      Var %d      ',i );
            for j = 1:p
                fprintf(1,'   %8.3f    ', PiHat(i,j) );
            end
            fprintf(1,'\n');
        end
        fprintf(1,  '-----------------------------------------------------------------------------------------------------\n');
        fprintf(1,  '                                                                                                     \n');
        fprintf(1,  '-----------------------------------------------------------------------------------------------------\n\n');
    end
end


if print2screen && printGammas && (k > 0)
    for l = 1:k
        GammaHatk = GammaHat( :, p*(l-1)+1 : p*l );
        GammaSEk = GammaSE( :, p*(l-1)+1 : p*l );
        
        fprintf(1,'    Lag Matrix %d (Gamma_%d):                                                                            \n', l, l );
        fprintf(1,  '-----------------------------------------------------------------------------------------------------\n');
        fprintf(1,    '      Variable  ' );
        for j = 1:p
            fprintf(1,    '       Var %d   ', j);
        end
        fprintf(1,'\n');
        fprintf(1,  '-----------------------------------------------------------------------------------------------------\n');
        for i = 1:p
            fprintf(1,    '      Var %d      ',i );
            for j = 1:p
                fprintf(1,'   %8.3f    ', GammaHatk(i,j) );
            end
            fprintf(1,'\n');
            fprintf(1,    '       SE %d       ',i );
            for j = 1:p
                fprintf(1,' (%8.3f  )  ', GammaSEk(i,j) );
            end
            fprintf(1,'\n');
        end
        fprintf(1,  '-----------------------------------------------------------------------------------------------------\n');
        fprintf(1,  'Note: Standard Errors in parenthesis.                                                                \n');
        fprintf(1,  '-----------------------------------------------------------------------------------------------------\n\n');
    end
end



% Plot characteristic roots to check for stationarity of process.
charPolyRoots = CharPolyRoots(dbHat(end), alphaHat, betaHat, GammaHat, print2screen)



% end
