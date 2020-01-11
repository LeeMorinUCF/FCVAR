function [ rankTestStats ] = RankTests(x, k, opt)
% function [ rankTestStats ] = RankTests(x, k, opt)
% Written by Michal Popiel and Morten Nielsen (This version 11.17.2014)
% Based on Lee Morin & Morten Nielsen (June 5, 2013)
%
% DESCRIPTION: Performs a sequence of  likelihood ratio tests
% 	for cointegrating rank.
%
% The results are printed to screen if the indicator print2screen is 1.
%
% input = vector or matrix x of data.
%       scalar k denoting lag length.
%       opt (object containing estimation options)
%
% output = rankTestStats structure with results from cointegrating rank
%           tests, containing the following (p+1) vectors with ith element
%           corresponding to rank = i-1:
%	dHat	(estimates of d)
%	bHat	(estimate of b)
%	LogL	(maximized log-likelihood)
%	LRstat  (LR trace statistic for testing rank r against rank p)
%	pv      (P-value of LR trace test, or "999" if P-value is not available)
%______________________________________________________


T = size(x, 1) - opt.N;
p = size(x, 2);

% Store user specified options for printing to screen because it will be
% turned off while looping over ranks.
tempPrint2Screen = opt.print2screen;

% Create output storage to be filled.
bHat   = zeros(p+1,1);
dHat   = zeros(p+1,1);
LogL   = zeros(p+1,1);
LRstat = zeros(p+1,1);
pv     = zeros(p+1,1);

% Do not print FCVAR estimation for each rank in the loop.
opt.print2screen = 0;

% Do not plot roots of characteristic polynomial for each lag in the loop.
opt.plotRoots = 0;

% Do not calculate standard errors.
opt.CalcSE = 0;

% For calculation of P-values
if(opt.rConstant || opt.levelParam)
    consT = 1;
else
    consT = 0;
end

% Estimate models for all ranks
for r = 0 : p
    [results] = FCVARestn(x,k,r, opt);
    dHat(r+1) = results.coeffs.db(1);
    bHat(r+1) = results.coeffs.db(2);
    LogL(r+1) = results.like;
end

% Calculate the LR statistics and P-values
for r = 0 : p-1

    LRstat(r+1) =  - 2*( LogL(r+1) - LogL(p+1) );

	p_val=[];
    % Get P-values, if
    % (1) no deterministic terms, or
    % (2) there is only restricted constant and d=b, or
    % (3) there is only a level parameter and d=b.
	if((~opt.rConstant && ~opt.unrConstant && ~opt.levelParam) ||...
       (opt.rConstant  && ~opt.unrConstant && opt.restrictDB) ||...
       (opt.levelParam && ~opt.unrConstant && opt.restrictDB) )
        p_val = get_pvalues(p-r, bHat(r+1), consT, LRstat(r+1), opt);
	end

    % If automatic calls to P-value program have not been installed or
    % enabled, then p_val is empty. Set it to 999 so that it can have a
    % value for storage in the rankTestStats matrix below.
    if(isempty(p_val))
        p_val = 999;
    end

    % Store P-values.
    pv(r+1) = p_val;
end

% Restore settings.
opt.print2screen = tempPrint2Screen;

% Print the results to screen.
if opt.print2screen

    % create a variable for output strings
    yesNo = {'No','Yes'};

    fprintf(1,'\n-----------------------------------------------------------------------------------------------------\n');
    fprintf(1,'                         Likelihood Ratio Tests for Cointegrating Rank                               \n');
    fprintf(1,'-----------------------------------------------------------------------------------------------------\n');
    fprintf(1,'Dimension of system:  %6.0f     Number of observations in sample:       %6.0f \n', p, T+opt.N);
    fprintf(1,'Number of lags:       %6.0f     Number of observations for estimation:  %6.0f \n', k, T);
    fprintf(1,'Restricted constant:  %6s     Initial values:                         %6.0f\n', yesNo{opt.rConstant+1}, opt.N );
    fprintf(1,'Unestricted constant: %6s     Level parameter:                        %6s\n', yesNo{opt.unrConstant+1}, yesNo{opt.levelParam+1} );
    fprintf(1,'-----------------------------------------------------------------------------------------------------\n');
    fprintf(1,'Rank \t  d  \t  b  \t Log-likelihood\t LR statistic\t P-value\n');
    for i = 1:p
        if (pv(i) ~= 999)
            fprintf(1,'%2.0f   \t%5.3f\t%5.3f\t%15.3f\t%13.3f\t%8.3f\n', i-1, dHat(i), bHat(i), LogL(i), LRstat(i), pv(i));
        else
            fprintf(1,'%2.0f   \t%5.3f\t%5.3f\t%15.3f\t%13.3f\t    ----\n', i-1, dHat(i), bHat(i), LogL(i), LRstat(i));
        end
    end
    fprintf(1,'%2.0f   \t%5.3f\t%5.3f\t%15.3f\t         ----\t    ----\n', i, dHat(i+1), bHat(i+1), LogL(i+1));
    fprintf(1,'-----------------------------------------------------------------------------------------------------\n');

end

% Return structure of rank test results.
rankTestStats.dHat   = dHat;
rankTestStats.bHat   = bHat;
rankTestStats.LogL   = LogL;
rankTestStats.LRstat = LRstat;
rankTestStats.pv     = pv;

end


function [pv] = get_pvalues(q, b, consT, testStat, opt)
% Written by Michal Popiel and Morten Nielsen (This version 10.22.2014)
%
% DESCRIPTION: This function calls the program FDPVAL in the terminal and
% 	returns the P-value based on the user's inputs. The function's
% 	arguments must be converted to strings in order to interact with the
% 	terminal.
%
% Input = q        (number of variables minus rank)
%         b        (parameter)
%         consT    (boolean variable indicating whether or not there is
%                   constant present)
%         testStat (value of the test statistic)
%	      opt (object containing estimation options)
% Output = pv (P-value for likelihood ratio test)
%_________________________________________________________________________

    if(b<0.5)
        % Series are stationary so use chi^2 with (p-r)^2 df, see JN2012
        pv = 1-chi2cdf(testStat, q^2);
    else
        % For non-stationary systems use simulated P-values from
		% MacKinnon and Nielsen (2014, Journal of Applied Econometrics)
		% and the C++ program conversion by Jason Rhinelander.

        % Convert user input to system commands and
        % translate arguments to strings
        args = sprintf('%g %g %g %g', q, b, consT, testStat);
        % Combine path, program, and arguments
        outCode = sprintf('%s %s', opt.progLoc, args);

        % Evaluate system command
        % Note: fdpval is a separately installed program.
        % For more information see: https://github.com/jagerman/fracdist
        % For download see https://github.com/jagerman/fracdist/releases
        [ flag , pval] = system([outCode]);

        % Note: string is returned, so it needs to be converted
        if(flag==0) % check if program was executed without errors
            pv = str2double(pval);
        else % program failed
            pv = [];
        end
    end
end
