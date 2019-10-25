function [ estimates ] = GetParams(x, k, r, db, opt)
% function [ estimates ] = GetParams(x, k, r, db, opt)
% Written by Michal Popiel and Morten Nielsen (This version 04.13.2016)
% Based on Lee Morin & Morten Nielsen (August 22, 2011)
% 
% DESCRIPTION: This function uses FWL and reduced rank regression to obtain
% the estimates of Alpha, Beta, Rho, Pi, Gamma, and Omega
%
% Input = x   (matrix of variables to be included in the system)
%         k   (number of lags)
%         r   (number of cointegrating vectors)
%         db  (value of d and b)
%         opt (object containing the estimation options)
% Output = estimates (Matlab structure containing the following)
%          - estimates.db (taken directly from the input)
%          - estimates.alphaHat
%          - estimates.betaHat
%          - estimates.rhoHat
%          - estimates.piHat
%          - estimates.OmegaHat
%          - estimates.GammaHat ( p x kp matrix [GammaHat1,...,GammaHatk])
%_________________________________________________________________________

    
    [ Z0, Z1, Z2, Z3 ] = TransformData(x, k, db, opt);
    
    T = size(Z0, 1);
    p = size(Z0, 2);
    p1 = size(Z1, 2);

    % --- Concentrate out the unrestricted constant if present --- %
    if(opt.unrConstant)
        % Note Z3*inv(Z3'*Z3)*Z3' is just a matrix of ones / T
        Z0hat = Z0 - Z3*( (Z3'*Z3)\Z3'*Z0 );
        Z1hat = Z1 - Z3*( (Z3'*Z3)\Z3'*Z1 );
        if(k>0)
            Z2hat = Z2 - Z3*( (Z3'*Z3)\Z3'*Z2 );
        else
            Z2hat = Z2;
        end
    else
        Z0hat = Z0;
        Z1hat = Z1;
        Z2hat = Z2;
    end
    
    % ---- FWL Regressions --- %    
    if (k == 0) 
        % No lags, so there are no effects of Z2.
        R0 = Z0hat;
        R1 = Z1hat;
    else
        % Lags included: Obtain the residuals from FWL regressions.
        R0 = Z0hat - Z2hat*( (Z2hat'*Z2hat)\Z2hat'*Z0hat );
        R1 = Z1hat - Z2hat*( (Z2hat'*Z2hat)\Z2hat'*Z1hat );
    end

       
    % Calculate Sij matrices for reduced rank regression.
    S00 = R0'*R0/T;
    S01 = R0'*R1/T;
    S10 = R1'*R0/T;
    S11 = R1'*R1/T;

    % Calculate reduced rank estimate of Pi.
    if r == 0
        betaHat = [];
        betaStar = [];
        alphaHat = [];
        PiHat = [];
        rhoHat = [];
        OmegaHat = S00;
  
    elseif ( r > 0 ) && ( r < p )
        [V,D] = eig( inv(S11)*S10*inv(S00)*S01 );
        V = sortrows( [ V' diag(D) ], p1+1 )';
        V = V(1:p1,:);
        betaStar = V( 1:p1, p1 : -1 : p1-r+1 );
       
        % If either alpha or beta is restricted, then the likelihood is
        %   maximized subject to the constraints. This section of the code
        %   replaces the call to the switching algorithm in the previous
        %   version.
        if(~isempty(opt.R_Alpha) || ~isempty(opt.R_Beta) )
            [ betaStar, alphaHat, OmegaHat ]...         
                = RstrctOptm_Switch(betaStar, S00, S01, S11, T, p, opt);
            betaHat = betaStar(1:p, 1:r);
            PiHat = alphaHat*betaHat';
        else
            % Otherwise, alpha and beta are unrestricted, but unidentified.
            alphaHat = S01*betaStar*inv(betaStar'*S11*betaStar);
            OmegaHat = S00 - alphaHat*betaStar'*S11*betaStar*alphaHat';
            betaHat = betaStar(1:p, 1:r);
            PiHat = alphaHat*betaHat';
            % Transform betaHat and alphaHat to identify beta.
            %   The G matrix is used to identify beta by filling the first 
            %   rxr block of betaHat with an identity matrix.
            G = inv(betaHat(1:r,1:r));
            betaHat = betaHat*G;
            betaStar = betaStar*G;
            % alphaHat is post multiplied by G^-1 so that Pi= a(G^-1)Gb' = ab'
            alphaHat = alphaHat*inv(G)';
        end
            

        % Extract coefficient vector for restricted constant model if required.
        if opt.rConstant
            rhoHat = betaStar( p1, : );
        else
            rhoHat = [];
        end
        
    else % (r = p) and do no need reduced rank regression
        V = S01*inv(S11);
        betaHat = V(:,1:p)';
        % For restriced constant, rho = last column of V.
        alphaHat = eye(p);
        PiHat = betaHat;
        OmegaHat = S00 - S01*inv(S11)*S10;

        % Extract coefficient vector for restricted constant model if required.
        if opt.rConstant
            rhoHat = V(:,p1)';
        else
            rhoHat = [];
        end
        betaStar = [ betaHat; rhoHat ];
    end

    % Calculate PiStar independently of how betaStar was estimated.
    PiStar = alphaHat*betaStar';

    % Perform OLS regressions to calculate unrestricted constant and
	%   GammaHat matrices if necessary.
    xiHat = [];
    if (k == 0)
        GammaHat = [];
    else
        if(r>0)
            GammaHat = ( inv(Z2hat'*Z2hat)*Z2hat'*(Z0hat - Z1hat*PiStar') )';
        else
            GammaHat = ( inv(Z2hat'*Z2hat)*Z2hat'*Z0hat )';
        end
    end

   if(opt.unrConstant==1)
       if(k>0)
           if(r>0)
               xiHat = ( inv(Z3'*Z3)*Z3'*(Z0 - Z1*PiStar' - Z2*GammaHat') )' ;
           else
               xiHat = ( inv(Z3'*Z3)*Z3'*(Z0 - Z2*GammaHat') )' ;
           end
       else
           if(r>0)
               xiHat = ( inv(Z3'*Z3)*Z3'*(Z0 - Z1*PiStar') )' ;
           else
               xiHat = ( inv(Z3'*Z3)*Z3'*Z0 )' ;
           end
       end
   end
   

    % --- Return the results in a structure --- %
    estimates.db = db;
    estimates.alphaHat = alphaHat;
    estimates.betaHat = betaHat;
    estimates.rhoHat = rhoHat;
    estimates.PiHat = PiHat;
    estimates.OmegaHat = OmegaHat;
    estimates.GammaHat = GammaHat;
    estimates.xiHat = xiHat;
    
end
