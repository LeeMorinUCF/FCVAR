function [betaStar, alphaHat, OmegaHat] = RstrctOptm_Switch(beta0, S00, S01, S11, T, p, opt)
% function [betaStar, alphaHat, OmegaHat] 
%                   = RstrctOptm_Switch(beta0, S00, S01, S11, T, p, opt)
% Written by Michal Popiel and Morten Nielsen (This version 03.29.2016)
% 
% DESCRIPTION: This function is imposes the switching algorithm of Boswijk
%   and Doornik (2004, page 455) to optimize over free parameters psi 
%   and phi directly, combined with the line search proposed by 
%	Doornik (2016, working paper). We translate between  (psi, phi) and 
%	(alpha, beta) using the relation of R_Alpha*vec(alpha) = 0 and 
%	A*psi = vec(alpha'), and R_Beta*vec(beta) = r_beta and 
%	H*phi+h = vec(beta). Note the transposes.
%
% Input = beta0 (unrestricted estimate of beta)
%         S00, S01, S11 (product moments)
%         T (number of observations)
%         p (number of variables)
%         opt (object containing the estimation options)
% Output = betaStar (estimate of betaStar)
%          alphaHat (estimate of alpha)
%          OmegaHat (estimate of Omega)
%_________________________________________________________________________

r  = size(beta0,2);
p1 = p+opt.rConstant;

% Restrictions on beta.
if(isempty(opt.R_Beta))
    H = eye(p1*r);
    h = zeros(p1*r,1);
else
    H = null(opt.R_Beta);
    h = opt.R_Beta'*inv(opt.R_Beta*opt.R_Beta')*opt.r_Beta;
end

% Restrictions on alpha.
%   We use the commutation matrix K_pr to transform vec(A) into vec(A'),
%   see Magnus & Neudecker (1988, p. 47, eqn (1)).
Ip  = eye(p);
Kpr = reshape(kron(Ip(:),eye(r)),p*r,p*r);
if(isempty(opt.R_Alpha))
    A = Kpr * eye(p*r); 
else
    A = null(opt.R_Alpha*inv(Kpr));
end

% Least squares estimator of Pi, used in calculations below.
PiLS = (S01*inv(S11))'; 
vecPiLS = PiLS(:);

% Starting values for switching algorithm.
betaStar = beta0;
alphaHat =  S01*beta0*inv(beta0'*S11*beta0);
OmegaHat = S00 - S01*betaStar*alphaHat' - alphaHat*betaStar'*S01' ...
        + alphaHat*betaStar'*S11*betaStar*alphaHat';

% Algorithm specifications.
iters   = opt.UncFminOptions.MaxFunEvals; 
Tol     = opt.UncFminOptions.TolFun;
conv    = 0;
i       = 0;

% Tolerance for entering the line search; epsilon_s in Doornik's paper.
TolSearch = 0.01;

% Line search parameters.
lambda = [1 1.2 2 4 8];
nS = length(lambda);
likeSearch = NaN(nS,1);
OmegaSearch = zeros(p,p,nS);

% Get candidate values for entering the switching algorithm.
vecPhi1 = inv(H'*kron(alphaHat'*inv(OmegaHat)*alphaHat,S11)*H)*...
        H'*(kron(alphaHat'*inv(OmegaHat),S11))*...
        (vecPiLS - kron(alphaHat, eye(p1))*h);

% Translate vecPhi to betaStar.
vecB = H*vecPhi1 + h;
betaStar = reshape(vecB,p1,r);    

% Candidate value of vecPsi.
vecPsi1 = inv(A'*kron(inv(OmegaHat),betaStar'*S11*betaStar)*A)*...
        A'*(kron(inv(OmegaHat),betaStar'*S11))*vecPiLS;

% Translate vecPsi to alphaHat.
vecA = A*vecPsi1; % This is vec(alpha')
alphaHat = reshape(inv(Kpr)*vecA,p,r);

% Candidate values of piHat and OmegaHat.
piHat1 = alphaHat*betaStar';
OmegaHat = S00 - S01*betaStar*alphaHat' - alphaHat*betaStar'*S01' ...
    + alphaHat*betaStar'*S11*betaStar*alphaHat';

% Calculate the likelihood.
like1 = - log(det(OmegaHat));

while(i<=iters && ~conv)
          
    % Update values for convergence criteria.
    piHat0  = piHat1;
    like0   = like1;

    if(i == 1)
        % Initialize candidate values. 
        vecPhi0_c = vecPhi1;
        vecPsi0_c = vecPsi1;
    end

    %  ---- alpha update step ---- %
    % Update vecPsi.
    vecPsi1 = inv(A'*kron(inv(OmegaHat),betaStar'*S11*betaStar)*A)*...
            A'*(kron(inv(OmegaHat),betaStar'*S11))*vecPiLS;
    
    % Translate vecPsi to alphaHat.
    vecA = A*vecPsi1; % This is vec(alpha')
    alphaHat = reshape(inv(Kpr)*vecA,p,r);    

    %  ---- omega update step ---- %
    % Update OmegaHat.
    OmegaHat = S00 - S01*betaStar*alphaHat' - alphaHat*betaStar'*S01' ...
        + alphaHat*betaStar'*S11*betaStar*alphaHat';
    
    %  ---- beta update step ---- %
    % Update vecPhi.
    vecPhi1 = inv(H'*kron(alphaHat'*inv(OmegaHat)*alphaHat,S11)*H)*...
            H'*(kron(alphaHat'*inv(OmegaHat),S11))*...
            (vecPiLS - kron(alphaHat, eye(p1))*h);
           
    % Translate vecPhi to betaStar.
    vecB = H*vecPhi1 + h;
    betaStar = reshape(vecB,p1,r);    
    
    %  ---- pi and likelihood update  ---- %
    % Update estimate of piHat.
    piHat1 = alphaHat*betaStar';
    
    % Update OmegaHat with new alpha and new beta.
    OmegaHat = S00 - S01*betaStar*alphaHat' - alphaHat*betaStar'*S01' ...
        + alphaHat*betaStar'*S11*betaStar*alphaHat';
    
    % Calculate the likelihood.
    like1 = - log(det(OmegaHat));

    if(i>0)
        
        % Calculate relative change in likelihood
        likeChange = (like1 - like0) / (1 + abs(like0));
        
        % Check relative change and enter line search if below tolerance.
        if(likeChange < TolSearch && opt.LineSearch)
            % Calculate changes in parameters.
            deltaPhi = vecPhi1 - vecPhi0_c;
            deltaPsi = vecPsi1 - vecPsi0_c;
            % Initialize parameter bins;
            vecPhi2 = NaN(length(vecPhi1), nS);
            vecPsi2 = NaN(length(vecPsi1), nS);
            % Values already calculated for lambda = 1.
            vecPhi2(:,1)       = vecPhi1; 
            vecPsi2(:,1)       = vecPsi1;
            likeSearch(1)      = like1;
            OmegaSearch(:,:,1) = OmegaHat;
            for iL = 2:nS
                % New candidates for parameters based on line search.
                vecPhi2(:,iL) = vecPhi0_c + lambda(iL)*deltaPhi;
                vecPsi2(:,iL) = vecPsi0_c + lambda(iL)*deltaPsi;
                
                % Translate to alpha and beta
                vecA = A*vecPsi2(:,iL); 
                alphaHat = reshape(inv(Kpr)*vecA,p,r);
                vecB = H*vecPhi2(:,iL) + h;
                betaStar = reshape(vecB,p1,r); 

                % Calculate and store OmegaHat
                OmegaSearch(:,:,iL) = S00 - S01*betaStar*alphaHat' - ...
                    alphaHat*betaStar'*S01' ...
                    + alphaHat*betaStar'*S11*betaStar*alphaHat';
                
                % Calculate and store log-likelihood
                likeSearch(iL) = - log(det(OmegaSearch(:,:,iL)));                
            end
            % Update max likelihood and OmegaHat based on line search.
            [ iSearch ] = find(likeSearch == max(max(likeSearch)));
            
            % If there are identical likelihoods, choose smallest
            % increment
            iSearch  = min(iSearch);
            like1    = likeSearch(iSearch);
            OmegaHat = OmegaSearch(:,:,iSearch);
            % Save old candidate parameter vectors for next iteration.
            vecPhi0_c  = vecPhi1;
            vecPsi0_c  = vecPsi1;
            
            % Update new candidate parameter vectors.
            vecPhi1  = vecPhi2(:,iSearch);
            vecPsi1  = vecPsi2(:,iSearch);
            % Update coefficients.
            vecA     = A*vecPsi1; 
            alphaHat = reshape(inv(Kpr)*vecA,p,r);
            vecB     = H*vecPhi1 + h;
            betaStar = reshape(vecB,p1,r); 
            % Update estimate of piHat.
            piHat1 = alphaHat*betaStar';
        end
               
        % Calculate relative change in likelihood
        likeChange = (like1 - like0) / (1 + abs(like0));
        
        % Calculate relative change in coefficients.
        piChange   = max(abs(piHat1(:) - piHat0(:)) ./ (1 + abs(piHat0(:))) );
               
        % Check convergence.        
        if(abs(likeChange) <= Tol && piChange <= sqrt(Tol))
            conv = 1;
        end
        
    end
    i=i+1;
end

end
