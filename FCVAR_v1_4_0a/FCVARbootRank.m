function [LRbs, H, mBS, mUNR] = FCVARbootRank(x, k, opt, r1, r2, B)
% function [LRbs, H, mBS, mUNR] = FCVARbootRank(x, k, opt, r1, r2, B)
% Written by Michal Popiel and Morten Nielsen (This version 08.06.2015)
% 
% DESCRIPTION: This function generates a distribution of a likelihood ratio
%           test statistic for the rank test using a Wild bootstrap, 
%			following the method of Cavaliere, Rahbek, and Taylor (2010). It 
%           takes the two ranks as inputs to estimate the model under the 
%           null and the model under the alternative.
%
% Input = x  (data - if k>0, actual data is used for initial values)
%         k  (number of lags)
%		  opt(estimation options)
%         r1 (rank under the null)
%         r2 (rank under the alternative)
%         B  (number of bootstrap samples)
% 
% Output = LRbs (B x 1 vector simulated likelihood ratio statistics)
%          pv (approximate p-value for LRstat based on bootstrap 
%                                                       distribution)
%          H (a Matlab structure containing LR test results, it is
%               identical to the output from HypoTest, with one addition,
%               namely H.pvBS which is the Bootstrap P-value)
%          mBS  (model estimates under the null)
%          mUNR (model estimates under the alternative)
%_________________________________________________________________________


% Calculate length of sample to generate, adjusting for initial values
T = size(x,1) - opt.N;

% Use first k+1 observations for initial values
data = x(1:k+1,:);

LR = zeros(B,1);

% Turn off output and calculation of standard errors for faster computation
opt.print2screen = 0;
opt.print2screen = 0;

mBS  = FCVARestn(x, k, r1, opt);
mUNR = FCVARestn(x, k, r2, opt);

H.LRstat = -2*(mBS.like - mUNR.like);

parfor j = 1:B
    
    % Display iteration count every 100 Bootstraps
    if(round((j+1)/10) == (j+1)/10)
        fprintf('iteration: %1.0f\n', j);
    end
    
    % (1) generate bootstrap DGP under the null
    xBS = FCVARsimBS(data, mBS, T);
    % append initial values to bootstrap sample
    BSs = [data; xBS];
    
    % (2) estimate unrestricted model
    mUNRbs =  FCVARestn(BSs, k, r2, opt);
    
    % (3) estimate restricted model (under the null)
    mRES =  FCVARestn(BSs, k, r1, opt);
    
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
    mUNR.like, mBS.like);
fprintf('Test results:\nLR statistic: \t %3.3f\nP-value (BS): \t %1.3f\n',...
    H.LRstat, H.pvBS);

end