function [LRbs, H, mBS, mUNR] = FCVARboot(x, k, r, optRES, optUNR, B)
% function [LRbs, H, mBS, mUNR] = FCVARboot(x, k, r, optRES, optUNR, B)
% Written by Michal Popiel and Morten Nielsen (This version 08.06.2015)
% 
% DESCRIPTION: This function generates a distribution of a likelihood ratio
%           test statistic using a Wild bootstrap, following the method of
%			Boswijk, Cavaliere, Rahbek, and Taylor (2013). It takes two sets 
%           of options as inputs to estimate the model under the null and the
%           unrestricted model. 
%
% Input = x      (data - if k>0, actual data is used for initial values)
%         k      (number of lags)
%         optRES (options object for restricted model under the null)
%         optUNR (options object to estimate unrestricted model)
%         B      (number of bootstrap samples)
% 
% Output = LRbs (B x 1 vector simulated likelihood ratio statistics)
%          pv   (approximate p-value for LRstat based on bootstrap 
%                                                       distribution)
%          H    (a Matlab structure containing LR test results, it is
%               identical to the output from HypoTest, with one addition,
%               namely H.pvBS which is the Bootstrap P-value)
%          mBS  (model estimates under the null)
%          mUNR (model estimates under the alternative)
%_________________________________________________________________________


% Calculate length of sample to generate, adjusting for initial values
T = size(x,1) - optRES.N;

% Use first k+1 observations for initial values
data = x(1:k+1,:);

LR = zeros(B,1);

% Turn off output and calculation of standard errors for faster computation
optUNR.print2screen = 0;
optRES.print2screen = 0;
optUNR.CalcSE = 0;
optRES.CalcSE = 0;

mBS  = FCVARestn(x, k, r, optRES);
mUNR = FCVARestn(x, k, r, optUNR);

fprintf('\nHypothesis test to bootstrap:\n');
H = HypoTest(mUNR, mBS);

% How often should the number of iterations be displayed
show_iters = 10;

parfor j = 1:B
    
    % Display iteration count every 100 Bootstraps
    if(round((j+1)/show_iters) == (j+1)/show_iters)
        fprintf('iteration: %1.0f\n', j);
    end
    
    % (1) generate bootstrap DGP under the null
    xBS = FCVARsimBS(data, mBS, T);
    % append initial values to bootstrap sample
    BSs = [data; xBS];
    
    % (2) estimate unrestricted model
    mUNRbs =  FCVARestn(BSs, k, r, optUNR);
    
    % (3) estimate restricted model (under the null)
    mRES =  FCVARestn(BSs, k, r, optRES);
    
    % (4) calculate test statistic   
    LR(j) = -2*(mRES.like - mUNRbs.like);

end

% Return sorted LR stats
LRbs = sort(LR);

% Calculate Bootstrap P-value (see ETM p.157 eq 4.62)
H.pvBS = sum(LRbs > H.LRstat)/B;

% Print output
fprintf('Bootstrap results:');
fprintf('\nUnrestricted log-likelihood: %3.3f\nRestricted log-likelihood:   %3.3f\n', ...
    H.loglikUNR, H.loglikR);
fprintf('Test results (df = %1.0f):\nLR statistic: \t %3.3f\nP-value: \t %1.3f\n',...
    H.df, H.LRstat, H.pv);
fprintf('P-value (BS): \t %1.3f\n', H.pvBS);

end