function [ betaHat, alphaHat, rhoHat, OmegaHat ]...
    = GetParamsSwitching(beta0, S00, S01, S11, T, p, r, EstnOptions)

% function [ betaHat, alphaHat, rhoHat, OmegaHat ]
%       = GetParamsSwitching(beta0, S00, S01, S11, T, p, r, EstnOptions)
% Lee Morin & Morten Nielsen
% August 21, 2011
%
% 
% GetParams(x, k, r, param, EstnOptions) returns the parameter estimates of
%       the fractional AR model using regression and reduced rank regression, 
%       by way of Johansen's switching algorithm.
% 
% input = matrix beta0 of starting values in the switching algorithm.
%       matrices S00, S01 and S11 of inner products, see Johansen (1996).
%       scalar T is the sample size.
%       scalar p is the dimension of the system.
%       scalar r is the cointegrating rank.
%       cell array EstnOptions.
% 
%   EstnOptions is an array of options for estimating the remaining parameters.
%   EstnOptions = { N, M, deterministics, R_Gamma, r_Gamma, [], [], [],
%                   nCols, A, H, a, b, h, [], [], [], [],
%                   [], [], [], print2screen, printGammas,
%                   restrictFminOptions }
%   The elements of EstnOptions are described in detail in the comments in
%       functions DefaultEstnOptions(db0, k, r, p),
%       and RestrictEstnOptions(db0, k, r, p, defaultFCVARoptions)
%       and AdjustEstnOptions(db0, k, r, p, initialFCVARoptions).
% 
% output = matrices of parameter estimates concentrated out of the likelihood function. 
%   estimates = { db, alphaHat, betaHat, rhoHat, PiHat, OmegaHat, GammaHat }.
%       GammaHat is reported as the p x kp matrix [ GammaHat1, ... , GammaHatk ].
%
% No function dependencies.
%______________________________________________________



% Determine desired set of options.
deterministics = EstnOptions(3);

% Determine the restrictions on the cointegrating relations.
% EstnOptions = { N, M, deterministics, R_Gamma, r_Gamma, R_CI, r_CI, CIparam0, 
%                   nCols, A, H, a, b, h, RestFminOptions };
nCols = EstnOptions{9};
nRestns = length(nCols);    % Note: nRestns = 0 when nCols = [].
A = EstnOptions{10};
H = EstnOptions{11};
a = EstnOptions{12};
b = EstnOptions{13};
h = EstnOptions{14};
Abar = EstnOptions{15};
Aperp = EstnOptions{16};
replaceCols = EstnOptions{17};
keepCols = EstnOptions{18};


% Determine the optimization settings for the switching algorithm.
RestFminOptions = EstnOptions{end};
MaxFunEvals = RestFminOptions.MaxFunEvals;
TolX = RestFminOptions.TolX;
TolFun = RestFminOptions.TolFun;



% Initialize betaHat.
betaHat = beta0;
for j = 1:nRestns
    if ~isempty(b(j)) % This(These) beta column(s) is (are) known.
        beta1 = b{j}; % Or: cell2mat(b(j))
    else
        % Estimate current beta (Note: Uses entire unrestricted beta0 matrix).
        [V,D] = eigs( beta0'*H{j}*inv(H{j}'*H{j})*H{j}'*beta0, beta0'*beta0 );
        vj = V(:, 1 : nCols(j) );
        beta1 = beta0*vj;
    end
    % Assemble updated beta0 matrix.
    betaHat(:, replaceCols{j} ) = beta1;
end

% Initialize alphaHat.
alphaHat = S01*betaHat*inv(betaHat'*S11*betaHat);
if ~isempty(a)
    for j = 1:nRestns
        if ~isempty(a(j))
            alphaHat(:, replaceCols{j} ) = a{j};
        end
    end
end



% Initialize stopping rules for switching algorithm.
i = 0;
% betaHat = beta0; % Already set in initialization loop.
betaHatLast = betaHat + 10*ones(size(betaHat));
logLike = 0;
logLikeLast = 10;

% Iterate estimation of beta and alpha.
while (  (i < MaxFunEvals) 
         && (norm(betaHat - betaHatLast) > TolX)...
         && (abs(logLike - logLikeLast) > TolFun)
      )
    
    % Update values used as stopping rules.
    i = i+1;
    betaHatLast = betaHat;
    logLikeLast = logLike;
    
    % Reset the inner product matrices.
    S00ab = S00;
    S01ab = S01;
    S11ab = S11;
    
    for j = 1:nRestns
        
        % Adjust Sij matrices for other columns known assumed known, if any.
        if nRestns > 1

            % Create matrix of other known betas, called beta0.
            beta0 = betaHat(:, keepCols{j} ); % Only used if nRestns > 1.
            % This is the matrix of other betas assumed known for now.
            % Same for alpha.
            alpha0 = alphaHat(:, keepCols{j});

            % I don't think I need these (alpha0 and beta0) otherwise.
            % Could call them knownCols.

            S00ab = S00ab - S01ab*beta0*alpha0' - alpha0*beta0'*S01ab' + alpha0*beta0'*S11ab*beta0*alpha0';
            S01ab = S01ab - alpha0*beta0'*S11ab;
        end
        
        
        % Correct Sij matrices for linear restriction on current alphaj.
        % If no restrictions, will have AjPerp = [] (Otherwise not defined).
        if ~isempty(AjPerp)
            % Linear restrictions on alpha (FWL out the extra regressors).
            S00ab = S00ab - S00ab*Aperp{j}*inv(Aperp{j}'*S00ab*Aperp{j})*Aperp{j}'*S00ab;
            S01ab = S01ab - S00ab*Aperp{j}*inv(Aperp{j}'*S00ab*Aperp{j})*Aperp{j}'*S01ab;
            S11ab = S11ab - S10ab*Aperp{j}*inv(Aperp{j}'*S00ab*Aperp{j})*Aperp{j}'*S01ab;
        end

        % Matrix of cross terms is symmetric.
        S10ab = S01ab'; 
        
        
        
        % Finally, it's time for estimation.
        if ( isempty(b(j)) && isempty(a(j)) ) 
            % Beta and alpha estimated without further restriction.
            [V,D] = eigs( H{j}'* S10ab * Abar{j} *inv(Abar{j}' *S00ab * Abar{j})* Abar{j}' *S01ab * H{j},...
                H{j}'* S11ab * H{j} );
            % Assemble updated beta0 matrix.
            phij = V(:,1:nCols(j));
            beta1 = H{j}*phij;
            betaHat(:, replaceCols{j} ) = H{j}*phij; % = beta1.
            
            alphaHat(:, replaceCols{j} ) = S01ab * beta1 * inv(beta1' * S11ab * beta1);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%.
            % Missing a minus sign?
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%.
            
        elseif ( ~isempty(b(j)) && isempty(a(j)) )  
            % This(These) beta column(s) is (are) known (but alpha unknown).
            beta1 = b{j};
            % Do nothing to betaHat since this known beta is already initialized correctly.
            
            alphaHat(:, replaceCols{j} ) = S01ab * beta1 * inv(beta1' * S11ab * beta1);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%.
            % Missing a minus sign?
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%.
            
            
        elseif ( isempty(b(j)) && ~isempty(a(j)) )  
            % This(These) alpha column(s) is (are) known (but beta unknown).
            
            % Generate a pseudo inverse of a(j) and apply to left-hand variable.
            aTildej = pinv(a(j)');
            
            betaHat(:, replaceCols{j} ) = aTildej * S01ab * inv(S11ab);
            
            
        end
        % Else both beta and alpha are known and are already set from initialization.
        
        
        
        
    end
    
    
    % Obtain estimate of Omega.
    OmegaHat = S00 - S01*betaHat*alphaHat' - alphaHat*betaHat'*S01' + alphaHat*betaHat'*S11*betaHat*alphaHat';
    
    % Calculate current value of likelihood.
    logLike = - T*p/2*( log(2*pi) + 1)  - T/2*log(det(OmegaHat));
    
    
end



% Define estimate of rhoHat, if necessary.
if strcmp(deterministics, 'restricted constant') 
    betaHat = betaHat(1:p,:);
    rhoHat = betaHat(p+1,:);
else
    rhoHat = [];
end



% end


