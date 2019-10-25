function [ estimates ] = GetParams(x, k, r, db, EstnOptions)
% 
% function [ estimates ] = GetParams(x, k, r, db, EstnOptions)
% Lee Morin & Morten Nielsen
% August 21, 2011
%
% 
% GetParams(x, k, r, db, EstnOptions) returns the parameter estimates of
%       the fractional AR model using regression and reduced rank regression.
% 
% input = vector or matrix of data x.
%       scalar k denotes the number of lags.
%       scalar r denotes the cointegrating rank.
%       vector or scalar db of fractional differencing parameters.
%           The size of db indicates whether or not restriction d = b is imposed
%           based on the number of free parameters.
%       cell array EstnOptions.
% 
%   EstnOptions is an array of options for estimating the remaining parameters.
%   EstnOptions = { N, M, deterministics, R_Gamma, r_Gamma, [], [], [],
%                   nCols, A, H, a, b, h, [], [], [], [],
%                   [], [], [], print2screen, printGammas,
%                   restrictFminOptions }
% The elements of EstnOptions are described in detail in the comments in
%       functions DefaultEstnOptions(db0, k, r, p),
%       and RestrictEstnOptions(db0, k, r, p, defaultFCVARoptions)
%       and AdjustEstnOptions(db0, k, r, p, initialFCVARoptions).
% 
% output = matrices of parameter estimates concentrated out of the likelihood function. 
%   estimates = { db, alphaHat, betaHat, rhoHat, PiHat, OmegaHat, GammaHat }.
%       GammaHat is reported as the p x kp matrix [ GammaHat1, ... , GammaHatk ].
% 
% Calls the functions TransformData(x, k, param, EstnOptions)
%   and GetParamsSwitching(betaHat, S00, S01, S11, T, p, r, EstnOptions)
%   if any restrictions are imposed on the estimation of alpha and beta.
% 
%______________________________________________________


% Determine desired set of options.
N = EstnOptions{1};
M = EstnOptions{2};
deterministics = EstnOptions(3);
R_Gamma = EstnOptions{4};
r_Gamma = EstnOptions{5};
nCols = EstnOptions{9};
nRestns = length(nCols);    % Note: nRestns = 0 when nCols = [].
A = EstnOptions{10};
H = EstnOptions{11};
a = EstnOptions{12};
b = EstnOptions{13};
h = EstnOptions{14};


% Transform data.
[ Z0, Z1, Z2 ] = TransformData(x, k, db, EstnOptions);
T = size(Z0, 1);
p = size(Z0, 2);

% Perform FWL regressions to concentrate out GammaHat matrices if necessary.
if k == 0
    R0 = Z0;
    R1 = Z1;
else
    Z2TZ2I = inv(Z2'*Z2);
    GammaHatFWL0 = Z2TZ2I*Z2'*Z0;
    GammaHatFWL1 = Z2TZ2I*Z2'*Z1;
    
    if (size(R_Gamma, 1) > 0) && (size(r_Gamma, 1) > 0)
        % Adjust GammaHatFWL matrices for restrictions.
        RestAdj = kron(eye(p), Z2TZ2I) * R_Gamma' * inv(R_Gamma*kron(eye(p), Z2TZ2I)*R_Gamma');
        
        vecGammaHatFWL0 = reshape(GammaHatFWL0, p*p*k, 1); 
        vecGammaHatFWL0 = vecGammaHatFWL0 + RestAdj*(r_Gamma - R_Gamma*vecGammaHatFWL0);
        GammaHatFWL0 = reshape(vecGammaHatFWL0, p*k, p);
        
        vecGammaHatFWL1 = reshape(GammaHatFWL1, p*p*k, 1); 
        vecGammaHatFWL1 = vecGammaHatFWL1 + RestAdj*(r_Gamma - R_Gamma*vecGammaHatFWL1);
        GammaHatFWL1 = reshape(vecGammaHatFWL1, p*k, p);
    end
    
    R0 = Z0 - Z2*GammaHatFWL0;
    R1 = Z1 - Z2*GammaHatFWL1;
end

% Calculate Sij matrices to calculate PiHat = alphaHat*betaHat'.
S00 = R0'*R0/T;
S01 = R0'*R1/T;
S10 = R1'*R0/T;
S11 = R1'*R1/T;

% Calculate reduced rank estimate of Pi.
if r == 0
    
    betaHat = 0;
    alphaHat = zeros(p,1);
    PiHat = zeros(p,p);
    OmegaHat = S00;
    
    % Extract coefficient vector for restricted constant model if required.
    if strcmp(deterministics, 'restricted constant') 
        rhoHat = 0;
    else
        rhoHat = [];
    end
    
    % For restricted constant, rho = 0 (defined from V later).
    % V = zeros(p+1, 1);
    
elseif ( r > 0 ) && ( r < p )
    
    [V,D] = eigs( inv(S11)*S10*inv(S00)*S01 );
    betaStar = V(:, 1:r);
    alphaHat = - S01*betaStar*inv(betaStar'*S11*betaStar);
    OmegaHat = S00 + alphaHat*betaStar'*S10;

    % Note: betaStar is the same regardless of the deterministics.
    betaHat = betaStar(1:p, 1:r);
    PiHat = alphaHat*betaHat';




    % Transform betaHat and alphaHat to identify beta.
    G = inv(betaHat(1:r,1:r));
    betaHat = betaHat*G;
    alphaHat = alphaHat*inv(G)';

    % Extract coefficient vector for restricted constant model if required.
    if strcmp(deterministics, 'restricted constant')
        rhoHat = V( p+1, 1: r );
        % Transform rhoHat to identify beta.
        rhoHat = rhoHat*G;
    else
        rhoHat = [];
    end
    
    

    % That was the initialization of the switching algorithm in case there are restrictions
    %   such as beta_i = H_i*phi_i or alpha_i = A_i*psi_i,
    %   or some columns of beta or alpha are known.
    if ( ~isempty(H) || ~isempty(A) || ~isempty(a) || ~isempty(b) || ~isempty(h) )
        [ betaHat, alphaHat, rhoHat, OmegaHat ]...
            = GetParamsSwitching(betaHat, S00, S01, S11, T, p, r, EstnOptions);
        PiHat = alphaHat*betaHat';
    end

    
else % (r = p)
    
    betaHat = eye(p);
    V = S01*inv(S11);
    % For restriced constant, rho = last column of V.
    alphaHat = V(:,1:p);
    PiHat = alphaHat;
    OmegaHat = S00 - S01*inv(S11)*S10;
    
    % Extract coefficient vector for restricted constant model if required.
    if strcmp(deterministics, 'restricted constant') 
        rhoHat = V(:,p+1)';
    else
        rhoHat = [];
    end
    
    % Note: Again, rhoHat could be split off betaStar for all cases at the end.
end




% Note: Could save space and if statement by defining PiStar with rho included.
%   Appending rho to PiHat will not change algorithm otherwise.

% Perform OLS regressions to calculate GammaHat matrices if necessary.
if strcmp(deterministics, 'restricted constant') 
    
    % dbLb1x = dbLb1x(:,1:p)*PiHat'  -  dbLb1x(:,p+1)*rhoHat*alphaHat';
    Z1 = Z1(:,1:p)*PiHat'  -  Z1(:,p+1)*rhoHat*alphaHat';
else
    % dbLb1x = dbLb1x(:,1:p)*PiHat';
    Z1 = Z1(:,1:p)*PiHat';
end

% Note: Might choose to concentrate and/or estimate Gamma in another program.

if k == 0
    GammaHat = [];
else
    % GammaHat = (xd - dbLb1x )'*dLbkx*inv(dLbkx'*dLbkx);
    % Are we missing a PiHat here?
    % GammaHat = (Z0 - Z1)'*Z2*inv(Z2'*Z2);
    
    % The transposed order of GammaHat seems more convenient for now.
    % Transpose back for display later.
    % Problem: Restrictions will correspond to Gamma in other order.
    
    % Z2TZ2I = inv(Z2'*Z2) defined above.
    GammaHat = Z2TZ2I*Z2'*(Z0 - Z1*PiHat');
    
    if (size(R_Gamma, 1) > 0) && (size(r_Gamma, 1) > 0)
        % Adjust GammaHat matrices for restrictions.
        % RestAdj = kron(eye(p), Z2TZ2I) * R_Gamma' * inv(R_Gamma*kron(eye(p), Z2TZ2I)*R_Gamma') defined above.
        
        vecGammaHat = reshape(GammaHat, p*p*k, 1); 
        vecGammaHat = vecGammaHat + RestAdj*(r_Gamma - R_Gamma*vecGammaHat);
        GammaHat = reshape(vecGammaHat, p*k, p);
        
    end
    
    GammaHat = GammaHat';
    
    
end



% Report parameter estimates as requested by user.
estimates = { db, alphaHat, betaHat, rhoHat, PiHat, OmegaHat, GammaHat };  






% end


