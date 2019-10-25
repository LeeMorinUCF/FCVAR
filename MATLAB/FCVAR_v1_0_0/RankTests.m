function [ rankTestStats ] = rankTests(x, k, db0, rankTestOptions)
% 
% function [ rankTestStats ] = rankTests(x, k, db0, rankTestOptions)
% Lee Morin & Morten Nielsen
% August 21, 2011
%
% rankTests(x, k, db0, options) performs a sequence of  likelihood ratio tests 
%       for cointegrating rank.
% 
% input = vector or matrix x of data.
%       scalar k denoting lag length.
%       vector or scalar db0 of starting values for d and b for the fractional differencing parameters.
%           The size of db0 indicates whether or not restriction d = b is imposed
%           based on the number of free parameters.
%       and cell array of options rankTestOptions output from DefaultEstnOptions(db0, k, r, p).
% 
% output = vector rankTestStats of cointegrating rank test statistics.
% 
% rankTestOptions = { dbEstnOptions, EstnOptions }, an array of 2 arrays of options 
%       for use in estimation routines.
% dbEstnOptions is an array of options for the estimation of the fractional differencing parameters d and b.
% dbEstnOptions = { constrained, C_db, c_db, R_db, r_db, dbFminOptions };
% EstnOptions is an array of options for estimating the remaining parameters.
%   EstnOptions = { N, M, deterministics, R_Gamma, r_Gamma, [], [], [],
%                   nCols, A, H, a, b, h, [], [], [], [],
%                   [], [], [], print2screen, printGammas,
%                   restrictFminOptions }
% The elements of EstnOptions are described in detail in the comments in
%       functions DefaultEstnOptions(db0, k, r, p)
%       and AdjustEstnOptions(db0, k, r, p, initialFCVARoptions).
% 
% 
% Calls the function FCVARlike(x, k, r, db, EstnOptions).
% 
%______________________________________________________


% Extract the cell arrays of options for estimation.
dbEstnOptions = rankTestOptions{1};
C_db = dbEstnOptions{2};
c_db = dbEstnOptions{3};
R_db = dbEstnOptions{4};
r_db = dbEstnOptions{5};
dbFminOptions = dbEstnOptions{end};
EstnOptions = rankTestOptions{2};
N = EstnOptions{1};
M = EstnOptions{2};
deterministics = EstnOptions{3};

% Indicator for printing results to screen.
print2screen = EstnOptions{22};


% Perform estimation sequence and calculate statistics.
T = size(x, 1) - N;
p = size(x, 2);
rankTestStats = zeros(p+1,5);

if print2screen
    fprintf(1,'\nPerforming estimation with full rank model.');
end

% Estimate full rank model first.
[ dbHat, fullRankLike, exitflag ] ...
        = fmincon(@( db ) -FCVARlike(x, k, p, db, EstnOptions), db0,...
        C_db, c_db, R_db, r_db, [], [], [], dbFminOptions );
% Record results for printing to screen.
% Expand dbHat into separate d and b if unconstrained.
dbHat = [ dbHat dbHat ];
dbHat = dbHat(1:2);
rankTestStats(p+1,:) = [ p dbHat fullRankLike 0 ];



if print2screen
    fprintf(1,'Successfully performed estimation with full rank model.\n');
end

% Estimate for all models of reduced rank.
for r = 0 : p-1
    
    
    
    if print2screen
        fprintf(1,'\nPerforming estimation with model of rank%d.\n', r);
    end
    
    [ dbHat, rRankLike, exitflag ] ...
        = fmincon(@( db ) -FCVARlike(x, k, r, db, EstnOptions), db0,...
        C_db, c_db, R_db, r_db, [], [], [], dbFminOptions );
        
        LRstat =  2*( rRankLike - fullRankLike );
        % Record results for printing to screen.
        % Expand dbHat into separate d and b if unconstrained.
        dbHat = [ dbHat dbHat ];
        dbHat = dbHat(1:2);
        rankTestStats(r+1,:) = [ r dbHat rRankLike LRstat ];
        
        if print2screen
            fprintf(1,'Successfully performed estimation with model of rank%d.\n', r);
        end
end


% Print the results to screen.

if print2screen

    fprintf(1,'\n\n-----------------------------------------------------------------------------------------------------\n');
    fprintf(1,'\n\n                         Likelihood Ratio Tests for Cointegrating Rank                               \n');
    fprintf(1,'\n\n-----------------------------------------------------------------------------------------------------\n\n');
    fprintf(1,'Dimension of System:  %8.0f     \t Number of observations for estimation: \t %8.0f \n', p, T);
    fprintf(1,'Deterministics:   %s \t Initial values:                          \t %8.0f\n', char(deterministics), N );
    fprintf(1,'Number of Lags:       %8.0f     \t Number of observations in sample:      \t %8.0f \n', k, T+N);
    if M > 0
        fprintf(1,'\n\nTruncation of fractional filter after %d lags.\n', M);
    end
    fprintf(1,'\n-----------------------------------------------------------------------------------------------------\n\n');
    fprintf(1,' \t Rank         d               b        Log-likelihood   LR statistic \n');
    fprintf(1,' \t  %2.0f \t %8.3f \t %8.3f \t %8.3f \t %8.3f \t \n', rankTestStats(1:p,:)');
    fprintf(1,' \t  %2.0f \t %8.3f \t %8.3f \t %8.3f \t    ----    \t \n', rankTestStats(p+1,1:4)');
    fprintf(1,'\n-----------------------------------------------------------------------------------------------------\n\n\n');
    
end



% end