function results = HypoTest(modelUNR, modelR)
% function results = HypoTest(modelUNR, modelR)
% Written by Michal Popiel and Morten Nielsen (This version 2.24.2015)
% 
% DESCRIPTION: This function performs a likelihood ratio test of the null
% 	hypothesis: "model is modelR" against the alternative hypothesis:
% 	"model is modelUNR".
%
% Input = modelUNR (structure of estimation results created for unrestricted model)
%         modelR (structure of estimation results created for restricted model)
% Output = results: a Matlab structure containing test results
%            - results.loglikUNR (loglikelihood of unrestricted model)
%            - results.loglikR   (loglikelihood of restricted model)
%            - results.df        (degrees of freedom for the test)
%            - results.LRstat    (likelihood ratio test statistic)
%            - results.p_LRtest  (P-value for test)
%_________________________________________________________________________

    % Calculate the test statistic.
    LR_test = 2*(modelUNR.like - modelR.like);

    % Calculate the degrees of freedom by taking the difference in free
    % parameters between the unrestricted and restricted model.
    df = modelUNR.fp - modelR.fp;
    
    % Calculate the P-value for the test.
    p_LRtest = 1-chi2cdf(LR_test, df);

    % Print output.
    fprintf('\nUnrestricted log-likelihood: %3.3f\nRestricted log-likelihood:   %3.3f\n', ...
        modelUNR.like, modelR.like);
    fprintf('Test results (df = %1.0f):\nLR statistic: \t %3.3f\nP-value: \t %1.3f\n',...
        df,LR_test,p_LRtest);

    % Return the test results in a Matlab structure.
    results.loglikUNR = modelUNR.like;
    results.loglikR   = modelR.like;
    results.df        = df;
    results.LRstat    = LR_test;
    results.pv        = p_LRtest;
    
end