function  LagSelect(x, kmax, r, order, opt )
% function LagSelect(x, kmax, r, order, opt )
% Written by Michal Popiel and Morten Nielsen (This version 3.31.2016)
% 
% DESCRIPTION: This program takes a matrix of variables and performs lag
% 	selection on it by using the likelihood ratio test. Output and test
% 	results are printed to the screen.
%
% Input = x     (matrix of variables to be included in the system)
%         kmax  (maximum number of lags)
%         r     (cointegration rank = number of cointegrating vectors)
%         order (order of serial correlation for white noise tests)
%         opt   (object containing estimation options)
% Output = none (only output to screen)
%_________________________________________________________________________


% Determine (initial) dimensions of system.
T = size(x, 1); % Length of sample (before truncation for initial values).
p = size(x, 2); % Number of variables.
printWN = 0;    % Set to zero to avoid printing output for each WN test

% Do not print FCVAR estimation for each lag in the loop.
opt.print2screen = 0;

% Do not plot roots of characteristic polynomial for each lag in the loop.
opt.plotRoots = 0;

% Do not calculate standard errors.
opt.CalcSE = 0;

% -------- Create storage bins for output ---------- %
% Row j of each of the following contains the associated results for
% 	lag length j+1
D        = zeros(kmax+1,2); % estimates of d and b
loglik   = zeros(kmax+1,1); % log-likelihood
LRtest   = zeros(kmax+1,1); % likelihood ratio test statistic for 
							%	significance of Gamma_{j+1}
pvLRtest = zeros(kmax+1,1); % likelihood ratio test P-value
aic      = zeros(kmax+1,1); % Akaike information criterion
bic      = zeros(kmax+1,1); % Bayesian information criterion
pvMVq    = zeros(kmax+1,1); % multivariate residual white noise Q-test P-value
pvWNQ    = zeros(kmax+1,p); % univariate residual white noise Q-test P-values
							%   for each residual
pvWNLM   = zeros(kmax+1,p); % univariate residual white noise LM-test P-values
							%   for each residual

  for k = 0:kmax
     
        % ----- Estimation ---------%         
        results = FCVARestn(x,k,r,opt);
       
        % ----- Record relevant output ---------%
        loglik(k+1) = results.like;
        D(k+1,:)    = results.coeffs.db;
        aic(k+1)    = -2*loglik(k+1)+ 2*results.fp;
        bic(k+1)    = -2*loglik(k+1)+ results.fp*log(T-opt.N);
              
        % ----- White noise tests ---------%
        [ ~, pvWNQ(k+1,:), ~, pvWNLM(k+1,:), ~, pvMVq(k+1,:) ]  = ...
                    mv_wntest(results.Residuals,order,printWN);   
		
		% ----- LR test of lag = k vs lag = k-1 -----%
		if k > 0
          LRtest(k+1)   = 2*(loglik(k+1)-loglik(k));
          pvLRtest(k+1) = 1-chi2cdf(LRtest(k+1), p^2);
		end
       
  end
  
  % Find lag corresponding to min of information criteria
  [~, i_aic] = min(aic);
  [~, i_bic] = min(bic);


  % --- Print output --- %
  % create a variable for output strings
  yesNo = {'No','Yes'};

  fprintf(1,'\n-----------------------------------------------------------------------------------------------------\n');
  fprintf(1,'                        Lag Selection Results \n');
  fprintf(1,'-----------------------------------------------------------------------------------------------------\n');
  fprintf(1,'Dimension of system:  %6.0f     Number of observations in sample:       %6.0f \n', p, T);
  fprintf(1,'Order for WN tests:   %6.0f     Number of observations for estimation:  %6.0f \n', order, T-opt.N);
  fprintf(1,'Restricted constant:  %6s     Initial values:                         %6.0f\n', yesNo{opt.rConstant+1}, opt.N );   
  fprintf(1,'Unrestricted constant: %6s     Level parameter:                        %6s\n', yesNo{opt.unrConstant+1}, yesNo{opt.levelParam+1} );
  fprintf(1,'-----------------------------------------------------------------------------------------------------\n');
  fprintf('k  r    d    b      LogL     LR    pv    AIC       BIC     pmvQ');
  for i=1:p
      fprintf(' pQ%1.0f  pLM%1.0f', i,i);
  end
   
  fprintf('\n');
  for k=kmax:-1:0
%       fprintf('%2.0f %2.0f %4.3f %4.3f %7.2f %6.2f %5.3f %8.2f %8.2f %4.2f',...          
%             k, r, D(k+1,:), loglik(k+1), LRtest(k+1),...
%           pvLRtest(k+1), aic(k+1), bic(k+1), pvMVq(k+1,:));
    fprintf('%2.0f %2.0f %4.3f %4.3f %7.2f %6.2f %5.3f %8.2f',...          
        k, r, D(k+1,:), loglik(k+1), LRtest(k+1),...
        pvLRtest(k+1), aic(k+1));
    % For AIC add asterisk if min value
    if(k+1 == i_aic); fprintf('*'); else fprintf(' '); end
    % Print BIC information criteria and add asterisk if min value
    fprintf(' %8.2f', bic(k+1));
    if(k+1 == i_bic); fprintf('*'); else fprintf(' '); end
    % Print multivariate white noise test P-values
    fprintf(' %4.2f', pvMVq(k+1,:));
    % Print the individual series white noise test P-values
    for i = 1:p
      fprintf(' %4.2f %4.2f', pvWNQ(k+1,i), pvWNLM(k+1,i));
    end
    fprintf('\n');
  end
  fprintf(1,'-----------------------------------------------------------------------------------------------------\n');
end

