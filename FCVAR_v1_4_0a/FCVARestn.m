function [ results ] = FCVARestn(x,k,r,opt)
% function [ results ] = FCVARestn(x,k,r,opt)
% Written by Michal Popiel and Morten Nielsen (This version 04.09.2016)
% 
% DESCRIPTION: This function performs estimation of the FCVAR system. It is
% 	the main function in the program with several nested functions, each
% 	described below. It estimates the model parameters, calculates the
% 	standard errors and the number of free parameters, obtains the residuals
% 	and the roots of the characteristic polynomial, and prints the output.
%
% Input = x (matrix of variables to be included in the system)
%         k (number of lags)
%         r (number of cointegrating vectors)
%         opt (object containing the estimation options)
% Output = results (a Matlab structure containing estimation results)
%            - results.startVals     (Starting values used for optimization)
%            - results.options       (Estimation options)
%            - results.like          (Model log-likelihood)
%            - results.coeffs        (Parameter estimates)
%            - results.rankJ         (Rank of Jacobian for 
%					identification condition)
%            - results.fp            (Number of free parameters)
%            - results.SE            (Standard errors)
%            - results.NegInvHessian (Negative of inverse Hessian matrix)
%            - results.Residuals     (Model residuals)
%            - results.cPolyRoots    (Roots of characteristic polynomial)
%_________________________________________________________________________

    global estimatesTEMP;
    
    % --- Preliminary steps --- %
    T = size(x,1) - opt.N; % number of observations
    p = size(x,2);         % number of variables
    
    % Update options based on initial user input.
    opt = updateRestrictions(opt,p,r);
    
    % Clear previous instances of coefficient estimates. This is done
    %   because estimatesTemp is a global structure that will not be cleared
    %   automatically if estimation is interrupted.
    estimatesTEMP = [];
      
	  
    % --- GRID SEARCH --- %
    
    % Hide all warnings for grid search.
    warning('off', 'all');
    
    % Perform grid search and store results as starting values for
    %   numerical optimization below.
    if(opt.gridSearch) 
        fprintf('\nRunning grid search over likelihood for k=%g, r=%g.\n',...
            k,r);
        fprintf('This computation can be slow.\n');
        fprintf('Set opt.gridSearch = 0 to skip it.\n');
        opt.db0 = LikeGrid(x,k,r,opt);
        % Change upper and lower bounds to limit the search to some small 
        %   interval around the starting values.
        opt.UB_db(1:2) = min(opt.db0(1:2) + [0.1 0.1], opt.dbMax);
        opt.LB_db(1:2) = max(opt.db0(1:2) - [0.1 0.1], opt.dbMin);
    else
        % Call to GetBounds returns upper/lower bounds for (d,b) or 
        %  depending on whether or not restrictions have been imposed. 
        [ opt.UB_db, opt.LB_db ] = GetBounds(opt);
        
        
    end
    
   
    % Turn warnings back on for main estimation.
    warning('on','all');
    
    % Hide warnings related to algorithm changes in numerical optimization.
    warning('off','optim:fminunc:SwitchingMethod');
    warning('off','optimL:fminunc:SwitchingMethod');
    
    % Hide warnings related to singular matrix in numerical optimization.
    warning('off','MATLAB:nearlySingularMatrix');

    
    % --- ESTIMATION --- %
	
    % Store equality restrictions, inequality restrictions, upper and lower
    %   bounds, and starting values for estimation because they need to be 
    %   adjusted based on the presence of level parameter or d=b restriction.
    Rpsi = opt.R_psi;
    rpsi = opt.r_psi;
    startVals = opt.db0(1:2); % starting values for mu get added later.
    Cdb = opt.C_db;
    cdb = opt.c_db;

    % If Rpsi is empty, then optimization is over (d,b), otherwise it is 
    %  over phi. If it is over phi, need to make adjustments to startVals 
    %  and make Cdb and cdb empty. 
    
    if(size(Rpsi,1)==1)
        H_psi = null(Rpsi);
        
        % Need to back out phi from given db0
        startVals = inv(H_psi'*H_psi)*H_psi'*startVals(1:2)';

        if(opt.gridSearch)
            % Translate from d,b to phi.
            UB = inv(H_psi'*H_psi)*H_psi'*opt.UB_db';
            LB = inv(H_psi'*H_psi)*H_psi'*opt.LB_db';      
        else
            % Otherwise GetBounds returns the values in terms of phi.
            UB = opt.UB_db;
            LB = opt.LB_db;
        end
            
        % Add warning about turning these off if non-empty?
        Cdb = [];
        cdb = [];
        
        % Turn off equality restrictions if only one restriction is imposed
        %  because they will be imposed inside the profile likelihood
        %  function that is being maximized.
        
        Rpsi = [];
        rpsi = [];
        
    else
        UB = opt.UB_db;
        LB = opt.LB_db;              
    end
       
    % If estimation involves level parameter, make appropriate 
    %   adjustments so that the dimensions of the restrictions correspond to
    %   the number of parameters being estimated, i.e. fill with zeros.
    if(opt.levelParam)
        % If there are equality restrictions add coefficients of 0 to the
        %  level parameter.
        if(~isempty(Rpsi))
            Rpsi = [Rpsi zeros(size(Rpsi,1),p)];      
        end
        % If there are inequality restrictions add coefficients of 0 to the
        %  level parameter.
        if(~isempty(Cdb))
            Cdb = [Cdb zeros(size(Cdb,1),p)];
        end
        % If grid search has not been used to set all initial values, add p
        %   starting values, set as the first observations on each variable,
        %   to the startVals variable to account for the level parameter.
        if(opt.gridSearch == 0)
            startVals = [startVals x(1,:)];
        else
            startVals = [startVals opt.db0(3:end)];
        end
        
        % Level parameter is unrestricted, but the length of UB and LB
        %   needs to match the number of parameters in optimization.
        UB = [UB  ones(1,p)*Inf];
        LB = [LB -ones(1,p)*Inf];
        
    end
    
        
   if(size(opt.R_psi)==2)
        % d,b are exactly identified by the linear restrictions and Rpsi is
        %  invertible. We use opt.R_psi here because Rpsi is adjusted
        %  depending on the presence of level parameters. Transpose is
        %  necessary to match input/output of other cases.
        dbTemp = (opt.R_psi'*inv(opt.R_psi*opt.R_psi')*opt.r_psi)';   
        y = x;
        if(opt.levelParam)
            % Optimize over level parameter for given (d,b).
            StartVal = y(1,:);
            [ muHat, maxLike, ~ ] ...
                = fminunc(@( params ) -FCVARlikeMu(x, dbTemp, params, k, r, opt), ...
                StartVal, opt.UncFminOptions );
            y = x - ones(T+opt.N,p)*diag(muHat);
        else
            maxLike = -FCVARlike(y, dbTemp, k, r, opt);
        end
        % Obtain concentrated parameter estimates.
        [ estimates ] = GetParams(y, k, r, dbTemp, opt);

        % Storing the estimates in a global structure here allows us to skip a
        %   call to GetParams after optimization to recover the coefficient
        %   estimates
        estimatesTEMP = estimates;
        % If level parameter is present, they are the last p parameters in the
        %   params vector
        if (opt.levelParam)
            estimatesTEMP.muHat = muHat;
        else
            estimatesTEMP.muHat = [];
        end
   else
       [ ~, maxLike, ~ ] ...
            = fmincon(@( params ) -FCVARlike(x, params, k, r, opt), ...
            startVals, Cdb, cdb, Rpsi, rpsi, LB, UB, [], opt.ConFminOptions );
   end
        
        
    % Store the updated estimation options.
    results.startVals = startVals;
    results.options = opt;
    results.options.UB_db = UB;
    results.options.LB_db = LB;        
        
        
    % Adjust the sign of the likelihood and store the results
    maxLike = -maxLike;    
    results.like = maxLike;
 
    % Coefficients are taken from a global defined in the likelihood
    %   function
    results.coeffs = estimatesTEMP;
 
	
    % ----- CHECK RANK CONDITION ---- %
    
    p1 = p + opt.rConstant;
    rankJ = []; %initialize the rank
    
    if(r>0) % If rank is zero then Alpha and Beta are empty
        
        if(isempty(opt.R_Beta))
            H_beta = eye(p1*r);
        else
            H_beta = null(opt.R_Beta);
        end
        
        % We use the commutation matrix K_pr to transform vec(A) into vec(A'),
		%   see Magnus & Neudecker (1988, p. 47, eqn (1)).
		Ip = eye(p);
        Kpr = reshape(kron(Ip(:),eye(r)),p*r,p*r);
        if(isempty(opt.R_Alpha))
            A = Kpr * eye(p*r);
        else
            A = null(opt.R_Alpha*inv(Kpr));
        end
		
        rA = size(A,2); % number of free parameters in alpha
        rH = size(H_beta,2); % number of free parameters in beta (including constant)

        % Following Boswijk & Doornik (2004, p.447) identification condition
        kronA = kron(eye(p), [results.coeffs.betaHat; results.coeffs.rhoHat])...
                                *[A zeros(p*r,rH)];
        kronH = kron(results.coeffs.alphaHat, eye(p1))*[zeros(p1*r,rA) H_beta];
        rankJ = rank(kronA + kronH);

        results.rankJ = rankJ;

    end
    
	
    % --- CHECK RANKS OF ALPHA AND BETA ---%
	
    % Check that alpha and beta have full rank to ensure that restrictions
    %   do not reduce their rank.
    if(rank(results.coeffs.alphaHat) < r)
       fprintf('\nWarning: Alpha hat has rank less than r!\n')
    end
    
    if( rank(results.coeffs.betaHat) < r)
        fprintf('\nWarning: Beta hat has rank less than r!\n')
    end
    
	
    % --- FREE PARAMETERS --- %
	
    % Compute the number of free parameters in addition to those in alpha
    %   and beta.
    [ fp ] = FreeParams(k, r, p, opt, rankJ);
    % Store the result.
    results.fp = fp;

	
    % --- STANDARD ERRORS --- %
    
    if(opt.CalcSE)
        % If any restrictions have been imposed, the Hessian matrix must be
        %   adjusted to account for them.
        if( ~isempty(opt.R_Alpha) ||~isempty(opt.R_psi) )

            % Create R matrix with all restrictions.

            % Count the number of restrictions on d,b. Note: opt.R_psi already 
            %  contains restrict DB, so the size() is only reliable if 
            %  it's turned off.
            if(~isempty(opt.R_psi))
                if(opt.restrictDB)
                    rowDB = size(opt.R_psi,1) - 1; 
                else
                    rowDB = size(opt.R_psi,1);
                end
            else
                % Otherwise d,b are unrestricted.
                rowDB = 0;
            end
            
            % Number of restrictions on alpha.
            rowA  = size(opt.R_Alpha,1);

            % Count the variables.
            colDB = 1 + ~opt.restrictDB;
            colA = p*r;
            colG = p*p*k;
            colMu = opt.levelParam*p;
            colRh = opt.unrConstant*p;
            % Length of vec(estimated coefficients in Hessian).
            R_cols  = colDB + colMu + colRh + colA + colG;

            % The restriction matrix will have rows equal to the number of
            %   restrictions.
            R_rows = rowDB + rowA;

            R = zeros(R_rows, R_cols);

            % Fill in the matrix R. 

            % Start with restrictions on (d,b) and note that if the model 
            %   with d=b is being estimated, only the first column of R_psi is
            %   considered. In that case, if there are zeros found in the first 
            %   column, the user is asked to rewrite the restriction. 
            if(rowDB>0)
                if(opt.restrictDB)
                    % If the model d=b is being estimated, only one
                    %  restriction can be imposed and that is on d.
                    R(1:rowDB,1:colDB) = opt.R_psi(1,1); 
                else                    
                    R(1:rowDB,1:colDB) = opt.R_psi; 
                end
            end
            % Put the R_Alpha matrix into the appropriate place in R.
            if(~isempty(opt.R_Alpha))
                R(1+rowDB: rowDB + rowA, ...
                    1 + colDB + colMu + colRh: colDB + colMu + colRh + colA)...
                    = opt.R_Alpha;
            end

            % Calculate unrestricted Hessian.
            [ H ] = FCVARhess(x, k, r, results.coeffs, opt);


            % Calculate the restricted Hessian.
            Q = -inv(H) + inv(H)*R'*inv(R*inv(H)*R')*R*inv(H);
        else
            % Model is unrestricted.
            [ H ] = FCVARhess(x, k, r, results.coeffs, opt);
            Q = -inv(H);
        end
    else
        NumCoeffs = length(SEmat2vecU(results.coeffs, k, r, p, opt)); 
        Q = zeros(NumCoeffs);
    end
    % Calculate the standard errors and store them.
    SE = sqrt(diag(Q));
    results.SE = SEvec2matU(SE,k,r,p, opt);
    results.NegInvHessian = Q;            

	
	
    % --- GET RESIDUALS --- %
    
    [ epsilon ] = GetResiduals(x, k, r, results.coeffs, opt);
    results.Residuals = epsilon;
    
	
    % --- OBTAIN ROOTS OF CHARACTERISTIC POLYNOMIAL --- %

    cPolyRoots = CharPolyRoots(results.coeffs, opt, k, r, p);
    results.cPolyRoots = cPolyRoots;
    
	
    % --- PRINT OUTPUT --- %
	
    if (opt.print2screen)
        if(~opt.CalcSE)
            fprintf('Warning: standard errors have not been calculated!\n');
        end
        
        % create a variable for output strings
        yesNo = {'No','Yes'};
        fprintf(1,'\n-----------------------------------------------------------------------------------------------------\n');
        fprintf(1,'                      Fractionally Cointegrated VAR: Estimation Results                              ');
        fprintf(1,'\n-----------------------------------------------------------------------------------------------------\n');
        fprintf(1,'Dimension of system:  %6.0f      Number of observations in sample:       %6.0f \n', p, T+opt.N);
        fprintf(1,'Number of lags:       %6.0f      Number of observations for estimation:  %6.0f \n', k, T);
        fprintf(1,'Restricted constant:  %6s      Initial values:                         %6.0f\n', yesNo{opt.rConstant+1}, opt.N );   
        fprintf(1,'Unrestricted constant:%6s      Level parameter:                        %6s\n', yesNo{opt.unrConstant+1}, yesNo{opt.levelParam+1} );
        if(size(opt.R_psi,1)==1)
            % 1 restriction.
            dbUB = H_psi*UB(1);
            dbLB = H_psi*LB(1);
            dbStart = H_psi*startVals(1);
            fprintf(1,'Starting value for d:    %1.3f    Parameter space for d: (%1.3f , %1.3f) \n', dbStart(1), dbLB(1), dbUB(1));
		    fprintf(1,'Starting value for b:    %1.3f    Parameter space for b: (%1.3f , %1.3f) \n', dbStart(2), dbLB(2), dbUB(2));
        else
            % Unrestricted or 2 restrictions.
            fprintf(1,'Starting value for d:    %1.3f    Parameter space for d: (%1.3f , %1.3f) \n', startVals(1), LB(1), UB(1));
		    fprintf(1,'Starting value for b:    %1.3f    Parameter space for b: (%1.3f , %1.3f) \n', startVals(2), LB(2), UB(2));
			fprintf(1,'Imposing d >= b:      %6s\n', yesNo{opt.constrained+1} );

        end
        fprintf(1,'-----------------------------------------------------------------------------------------------------\n');
        fprintf(1,'Cointegrating rank:   %10.0f  AIC:            %10.3f \n', r, -2*maxLike + 2*fp);
        fprintf(1,'Log-likelihood:       %10.3f  BIC:            %10.3f \n', maxLike, -2*maxLike + fp*log(T));
        fprintf(1,'log(det(Omega_hat)):  %10.3f  Free parameters:%10.0f \n', log(det(results.coeffs.OmegaHat)), fp);
        fprintf(1,'-----------------------------------------------------------------------------------------------------\n');
        fprintf(1,    '    Fractional parameters:                                                                             \n');
        fprintf(1,    '-----------------------------------------------------------------------------------------------------\n');
        fprintf(1,    '    Coefficient              \t Estimate              \t  Standard error \n');
        fprintf(1,    '-----------------------------------------------------------------------------------------------------\n');
        fprintf(1,    '         d                   \t %8.3f              \t     %8.3f                \n', results.coeffs.db(1), results.SE.db(1));
        if ~opt.restrictDB
            fprintf(1,'         b                   \t %8.3f              \t     %8.3f                \n', results.coeffs.db(2), results.SE.db(2));
        end
        fprintf(1,  '-----------------------------------------------------------------------------------------------------\n');
        fprintf(1,  '-----------------------------------------------------------------------------------------------------\n');
        if r > 0
            if opt.rConstant
                varList = '(beta and rho):';
            else
                varList = '(beta):        ';
            end
            fprintf(1,'    Cointegrating equations %s                                                          \n', varList);
            fprintf(1,  '-----------------------------------------------------------------------------------------------------\n');
            fprintf(1,    '      Variable      ' );
            for j = 1:r
                fprintf(1,    '  CI equation %d  ', j);
            end
            fprintf(1,'\n');
            fprintf(1,  '-----------------------------------------------------------------------------------------------------\n');
            for i = 1:p
                fprintf(1,    '        Var%d       ',i );
                for j = 1:r
                    fprintf(1,'    %8.3f     ', results.coeffs.betaHat(i,j) );
                end
                fprintf(1,'\n');
            end
            if opt.rConstant
                fprintf(1,    '      Constant     ' );
                for j = 1:r
                    fprintf(1,'    %8.3f     ', results.coeffs.rhoHat(j) );
                end
                fprintf(1,'\n');
            end
            fprintf(1,  '-----------------------------------------------------------------------------------------------------\n');
            if (isempty(opt.R_Alpha) && isempty(opt.R_Beta) )
                fprintf(1,  'Note: Identifying restriction imposed.                                                               \n');
                fprintf(1,  '-----------------------------------------------------------------------------------------------------\n');
            end
            fprintf(1,'    Adjustment matrix (alpha):                                                                         \n' );
            fprintf(1,  '-----------------------------------------------------------------------------------------------------\n');
            fprintf(1,    '      Variable      ' );
            for j = 1:r
                fprintf(1,    '  CI equation %d  ', j);
            end
            fprintf(1,'\n');
            fprintf(1,  '-----------------------------------------------------------------------------------------------------\n');
            for i = 1:p
                fprintf(1,    '        Var %d      ',i );
                for j = 1:r
                    fprintf(1,'    %8.3f     ', results.coeffs.alphaHat(i,j) );
                end
                fprintf(1,'\n');
                fprintf(1,    '         SE %d      ',i );
                    for j = 1:r
                        fprintf(1,'   (%8.3f  )  ', results.SE.alphaHat(i,j) );
                    end
                fprintf(1,'\n');
            end
            fprintf(1,  '-----------------------------------------------------------------------------------------------------\n');
            fprintf(1,  'Note: Standard errors in parenthesis.                                                                \n');
            fprintf(1,  '-----------------------------------------------------------------------------------------------------\n');
            fprintf(1,'    Long-run matrix (Pi):                                                                       \n' );
            fprintf(1,  '-----------------------------------------------------------------------------------------------------\n');
            fprintf(1,    '      Variable  ' );
            for j = 1:p
                fprintf(1,    '       Var %d   ', j);
            end
            fprintf(1,'\n');
            fprintf(1,  '-----------------------------------------------------------------------------------------------------\n');
            for i = 1:p
                fprintf(1,    '      Var %d      ',i );
                for j = 1:p
                    fprintf(1,'   %8.3f    ', results.coeffs.PiHat(i,j) );
                end
                fprintf(1,'\n');
            end
            fprintf(1,  '-----------------------------------------------------------------------------------------------------\n');
        end
    end

% Print level parameter if present.
if (opt.print2screen && opt.levelParam)
    fprintf(1,  '\n-----------------------------------------------------------------------------------------------------\n');
    fprintf(1,'    Level parameter (mu):                                                                         \n' );
    fprintf(1,  '-----------------------------------------------------------------------------------------------------\n');
    for i = 1:p
        fprintf(1,    '        Var %d      ',i );
        fprintf(1,'    %8.3f     ', results.coeffs.muHat(i) );
        fprintf(1,'\n');
        fprintf(1,    '         SE %d      ',i );
        fprintf(1,'   (%8.3f  )  ', results.SE.muHat(i) );
        fprintf(1,'\n');
    end
    fprintf(1,  '-----------------------------------------------------------------------------------------------------\n');
    fprintf(1,  'Note: Standard errors in parenthesis (from numerical Hessian) but asymptotic distribution is unknown. \n');
    fprintf(1,  '-----------------------------------------------------------------------------------------------------\n');
end

% Print unrestricted constant if present.
if (opt.print2screen && opt.unrConstant)
    fprintf(1,  '\n-----------------------------------------------------------------------------------------------------\n');
    fprintf(1,'    Unrestricted constant term:                                                                     \n' );
    fprintf(1,  '-----------------------------------------------------------------------------------------------------\n');
    for i = 1:p
        fprintf(1,    '        Var %d      ',i );
        fprintf(1,'    %8.3f     ', results.coeffs.xiHat(i) );
        fprintf(1,'\n');
        fprintf(1,    '         SE %d      ',i );
        fprintf(1,'   (%8.3f  )  ', results.SE.xiHat(i) );
        fprintf(1,'\n');
    end
    fprintf(1,  '-----------------------------------------------------------------------------------------------------\n');
    fprintf(1,  'Note: Standard errors in parenthesis (from numerical Hessian) but asymptotic distribution is unknown. \n');
    fprintf(1,  '-----------------------------------------------------------------------------------------------------\n');
end    
    
% Print Gamma coefficients if required.
if opt.print2screen && opt.printGammas && (k > 0)
    for l = 1:k
        GammaHatk = results.coeffs.GammaHat( :, p*(l-1)+1 : p*l );
        GammaSEk = results.SE.GammaHat( :, p*(l-1)+1 : p*l );
 
        fprintf(1,'    Lag matrix %d (Gamma_%d):                                                                            \n', l, l );
        fprintf(1,  '-----------------------------------------------------------------------------------------------------\n');
        fprintf(1,    '      Variable  ' );
        for j = 1:p
            fprintf(1,    '       Var %d   ', j);
        end
        fprintf(1,'\n');
        fprintf(1,  '-----------------------------------------------------------------------------------------------------\n');
        for i = 1:p
            fprintf(1,    '      Var %d      ',i );
            for j = 1:p
                fprintf(1,'   %8.3f    ', GammaHatk(i,j) );
            end
            fprintf(1,'\n');
            fprintf(1,    '       SE %d       ',i );
            for j = 1:p
                fprintf(1,' (%8.3f  )  ', GammaSEk(i,j) );
            end
            fprintf(1,'\n');
        end
        fprintf(1,  '-----------------------------------------------------------------------------------------------------\n');
        fprintf(1,  'Note: Standard errors in parenthesis.                                                                \n');
        fprintf(1,  '-----------------------------------------------------------------------------------------------------\n');
    end
end

% Print roots of characteristic polynomial if required.
if (opt.print2screen && opt.printRoots)
    fprintf(1,  '\n-----------------------------------------------------------------------------------------------------\n');
    fprintf(1,  '    Roots of the characteristic polynomial                                                           \n');
    fprintf(1,  '-----------------------------------------------------------------------------------------------------\n');
    fprintf(1,  '    Number     Real part    Imaginary part       Modulus                                             \n');
    fprintf(1,  '-----------------------------------------------------------------------------------------------------\n');
    for j = 1:length(cPolyRoots)
        fprintf(1, '      %2.0f       %8.3f       %8.3f         %8.3f                                        \n',...
            j, real(cPolyRoots(j)), imag(cPolyRoots(j)), abs(cPolyRoots(j)) );
    end
    fprintf(1,  '-----------------------------------------------------------------------------------------------------\n\n');
end

% Print notifications regarding restrictions.
if (opt.print2screen &&  (~isempty(opt.R_Alpha) || ~isempty(opt.R_psi) ...
        || ~isempty(opt.R_Beta) ))
    fprintf(1,  '\n-----------------------------------------------------------------------------------------------------\n');
    fprintf(1,  'Restrictions imposed on the following parameters:\n');
    if(~isempty(opt.R_psi))
        fprintf('- Psi. For details see "options.R_psi"\n');
    end
    if(~isempty(opt.R_Alpha))
        fprintf('- Alpha. For details see "options.R_Alpha"\n');
    end
    if(~isempty(opt.R_Beta))
        fprintf('- Beta. For details see "options.R_Beta"\n');
    end
    fprintf(1,  '-----------------------------------------------------------------------------------------------------\n\n');
end

end










