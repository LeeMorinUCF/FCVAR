classdef EstOptions  
% classdef EstOptions 
% Written by Michal Popiel and Morten Nielsen (This version 04.11.2016)
% 
% DESCRIPTION: This class defines the estimation options used in the FCVAR
% 	estimation procedure and the related programs. Assigning this class
% 	to a variable stores the default properties defined below in that
% 	variable. In addition to the properties, the methods section includes the
% 	function updateRestrictions which performs several checks on the
% 	user-specified options prior to estimation.
%______________________________________________________

    properties
        
        % --- ESTIMATION OPTIONS --- %
        % Estimation options for unconstrained optimization.
        UncFminOptions = optimset('MaxFunEvals', 1000, 'TolFun', 1e-8,...
                        'TolX', 1e-8, 'Display', 'off');
        
        % Estimation options for constrained optimization.
		ConFminOptions = optimset('MaxFunEvals', 1000, 'TolFun', 1e-8,...
                        'TolX', 1e-8, 'Display', 'off', ...
                        'Algorithm', 'interior-point');
        
        % Activate live search for switching algorithm in restricted model
        % estimation.
        LineSearch = 1;
        
        % Variation on the grid search to find local max. 
        LocalMax   = 1;
                          
        % Set upper and lower bound for d,b parameters.
        dbMax = 2;
        dbMin = 0.01;

        % Starting value for optimization.
        db0 = [1 1]; 
        
        % Constrain parameter space for d,b. If set to 1, then impose
        % dbMin<= b <= d <= dbMax. If this option is set to 0 then 
        % dbMin<= d <= dbMax and dbMin <= b <= dbMax are imposed separately.
        constrained = 1;
        
        % If restrictDB = 1, then d=b is imposed. If = 0, the restriction 
		% is not imposed.
        restrictDB = 1;
        
        % Number of observations reserved for initial values to be
        % conditioned upon in the estimation.
        N = 0;     
        
		
        % --- DETERMINISTICS --- %
        % Unrestricted constant.
        unrConstant = 0;
        
        % Restricted constant.
        rConstant   = 0;
        
        % Level parameter.
        levelParam  = 1;
        
		
        % --- RESTRICTIONS --- %
        % Inequality constraints on parameters d and b.
        % Specified as C_db * [d; b] <= c_db
        % Note: these restrictions are non-standard and should only be
        % specified if constrained=0. Furthermore, the grid search is not
        % equipped to handle these types of restrictions.
        C_db = [];
        c_db = [];
        
        % Upper and lower bounds for parameters d and b.
        % Note: these options are set automatically in the
        % updateRestrictions function below based on inputs of 'dbMax',
        % and 'dbMin'.
        UB_db = [];
        LB_db = [];
        
        % Equality constraints on parameters d and b.
        % Specified as R_psi * [d; b] = r_psi
        R_psi = [];
        r_psi = [];
        
              
        % Restrictions on Alpha matrix. 
        % Specified as R_Alpha*vec(Alpha) = r_Alpha
        % Note: r_Alpha can only have 0's
        R_Alpha = [];
        r_Alpha = []; 
        
        % Restrictions on Beta matrix.
        % Specified as R_Beta*vec(Beta) = r_Beta
        % Note: r_Beta can have non-zero elements.
        R_Beta = [];
        r_Beta = [];


        % --- OUTPUT OPTIONS --- %
        % If set to 0, print2screen prevents all output from being printed.
        print2screen = 1;
        
        % If set to 0, printGammas prevents the estimated coefficients in
        % the Gamma matrix from being printed.
        printGammas = 1;
        
        % If set to 0, printRoots prevents the roots of the
        % characteristic polynomial from being printed.
        printRoots = 1;
        
        % If set to 0, plotRoots prevents the roots of the
        % characteristic polynomial from being plotted.
        plotRoots = 1;

		
        % --- Grid search --- %
        % gridSearch set to 1 makes the program evaluate the likelihood
        % over a grid of (d,b) and choose the max from the grid as the 
		% starting value for the numerical optimization. It is computationally
		% costly to perform but sometimes yields more accurate results.
        gridSearch = 1;
        
        % plotLike makes the program plot the likelihood function over the
        % grid of (d,b) when the grid search is selected.
        plotLike = 1;
        
        % progress prints the progress of the grid search in
        % either a waitbar window (=1) or to the commandline (=2) or not at
        % all (=0).
        progress = 1;
		
        % If progress ~= 0, print progress every updateTime seconds.
        updateTime = 5;
        
        % --- Set location path with program name --- %
		% Location path with program name of fracdist program, if installed,
		%   for calculation of P-values, see MacKinnon and Nielsen (2014, JAE).
		% If this is not installed, the option can be left as-is or blank.
        % Note: use both single (outside) and double quotes (inside). This
        %   is especially important if path name has spaces.
		
        % Linux example:
		progLoc = '"/usr/bin/fdpval"';  
		
        % Windows example: 
		%   program located in folder fdpval in current directory
		% progLoc = '".\fdpval\fdpval"';    
        
        % Calculate SE's. This option should not be changed by the user.
        %   Calculating standard errors is omitted in the LagSelection.m
        %   RankTests.m functions to speed up estimation.
        CalcSE = 1;

    end
    
    methods
        function newObj = updateRestrictions(obj, p, r)
            % function newObj = updateRestrictions(obj, p, r)
            % Description: This function is called prior to estimation and
            % 	performs several checks to verify the input in the estimation
            % 	options:
            % 	- Check that only one option for deterministics is chosen,
            % 	- Check that starting values have been specified correctly,
            % 	- Update dbMin, dbMax, based on user-specified options,
            % 	- Check for appropriate dimensions and redundancies in the
            % 		restriction matrices R_psi, R_Alpha and R_Beta.
            % Input = obj (estimation option class variable)
            %         p (number of variables in the system)
            %         r (cointegrating rank)
            % Output = newObj (updated object)
            % Dependencies: called by FCVARestn().

            % --- Deterministics --- %
            
            % Check for redundant deterministic specifications.
            if(obj.levelParam && obj.rConstant)
                obj.rConstant = 0;
                fprintf('\nWarning, both restricted constant and level parameter selected.\n');
                fprintf('Only level parameter will be included.\n');
            end


            % --- Fractional Parameters d,b --- %
			
            % Adjust starting values. 
            
            % If d=b is imposed but the starting values are different, then
            % adjust them so that they are the same for faster and more
            % accurate computation.
            if(obj.restrictDB && length(obj.db0)>1 &&...
                    obj.db0(1)~=obj.db0(2))
                obj.db0 = [obj.db0(1) obj.db0(1)];
                fprintf('\nWarning, d=b imposed but starting values are different for d,b.\n');
                fprintf('db0 set to [%g %g].\n', obj.db0);
            end
            
            % if only one starting value is specified, then d0 = b0. 
            if(length(obj.db0)<2)
                obj.db0 = [obj.db0 obj.db0];
            end

            % Check if too many parameters specified in dbMin/dbMax
            if ( length(obj.dbMin) > 2 || length(obj.dbMax) > 2 )
                 fprintf('\nWARNING: Too many parameters specified in dbMin or dbMax.\n');
                 fprintf('\nOnly the first two elements will be used to set bounds in estimation.\n');
                 % Cut off extra bounds
                 if ( length(obj.dbMin) > 2 )
                     obj.dbMin = obj.dbMin(1:2);
                 end
                 if ( length(obj.dbMax) > 2 )
                     obj.dbMax = obj.dbMax(1:2);
                 end
            end

            % Assign the same min/max values for d,b if only one set of
            % bounds is provided. This is mostly for backwards
            % compatibility with previous versions. 
            % Minimum
            if( length(obj.dbMin) == 1 )               
                obj.dbMin = [obj.dbMin obj.dbMin];
                % Notification to user.
                fprintf('\nWarning: minimum value specified for d only.\n');
                fprintf('Min value of b set to %2.4f\n', obj.dbMin(2));
            else
                % Set as row vector, regardless of input.
                obj.dbMin = reshape(obj.dbMin,1,2);
            end
            % Maximum
            if(length(obj.dbMax) == 1)
                obj.dbMax = [obj.dbMax obj.dbMax];
                % Notification to user.
                fprintf('\nWarning: maximum value specified for d only.\n');
                fprintf('Max value of b set to %2.4f\n', obj.dbMax(2));
            else
                % Set as row vector, regardless of input.               
                obj.dbMax = reshape(obj.dbMax,1,2);
            end
                
            % Check for consistency among min and max values for fractional
            % parameters.
            if( obj.dbMin(1) > obj.dbMax(1) )
                fprintf('\nWarning, min > max for bounds on d.\n');
                error('Invalid bounds inmposed on d.');
            end            
            if( obj.dbMin(2) > obj.dbMax(2) )
                fprintf('\nWarning, min > max for bounds on b.\n');
                error('Invalid bounds inmposed on b.');
            end          
            
            % If both inequality constraints and d>=b is imposed then
            % inconsistencies could arise in the optimization. Furthermore,
            % grid search is not equipped to deal with other types of
            % inequalities. 
            if(~isempty(obj.C_db) && obj.constrained)
                fprintf('\nWarning, inequality restrictions imposed and constrained selected.\n');
                fprintf('Restrictions in C_db override constrained options.\n');
                obj.constrained = 0;
                if(obj.gridSearch)
                    fprintf('Grid search has been switched off.\n');
                    obj.gridSearch = 0;
                end
            end
           
            % If restrictDB is selected, then only restrictions on d
            % are allowed. 
            if (~isempty(obj.R_psi))
                % First check if restriction is valid
                if(obj.restrictDB && ((obj.R_psi(1,1) == 0) || ...
                        (obj.R_psi(1,2) ~= 0)) )
                    fprintf('\nError in R_psi. When restrictDB = 1, only ');
                    fprintf('restrictions on d can be imposed.\n');
                    error('Invalid restriction for d=b model.');                   
                end
                % If non-zero restrictions haven't been specified by
                % the user, then set r_psi to zero's.
                if(isempty(obj.r_psi))
                    obj.r_psi = 0;
                end
            end
            
            % If d=b is imposed, add it to the other equality restrictions on d,b.
            if(obj.restrictDB)
                obj.R_psi = [obj.R_psi; [1 -1]];
                obj.r_psi = [obj.r_psi; 0];
                if(obj.constrained)
                    fprintf('\nNote: Redundant options. Both constrained (d>=b) and restrict (d=b) selected.');
                    fprintf('\n Only d=b imposed.\n');
                    % Turn off (d>=b) constraint.
                    obj.constrained = 0;
                end
            end
            
            % Check restrictions on fractional parameters.
            if (isempty(obj.R_psi))
	            % Ensure that parameter space for fractional parameters is non-empty.
                if(obj.dbMax(1) < obj.dbMin(2) && obj.constrained)
                    fprintf('\nWarning: Redefine restrictions on fractional parameters.\n')
                    error('Empty parameter space for (d,b).')
	            end
			else
               % Ensure that parameter space for fractional parameters is non-empty.
                [ UB, LB ] = GetBounds(obj);
                if(LB > UB)
                    fprintf('\nWarning: Redefine restrictions on fractional parameters.\n')
                    error('Empty parameter space for (d,b).')
                end
                % Check if grid search is necessary.
                if ( (size(obj.R_psi,1)>1 && obj.gridSearch)  )
                    fprintf('\nd and b are exactly identified by imposed restrictions\n');
                    fprintf('so grid search has been turned off.\n');
                    obj.gridSearch = 0;			
                end            
                % Check for redundancies.
                if(rank(obj.R_psi) < size(obj.R_psi,1))
                    fprintf('\nWARNING: R_psi has reduced rank!\n');
                    fprintf('\nRedefine R_psi with linearly independent restrictions.\n');
                    error('Redundant restrictions in R_psi matrix');                  
                end   
            end                        
            
			            
            % --- Alpha and Beta --- %
            
            % Define p1 to be number of rows of betaStar.
            p1 = p + obj.rConstant;
            
            % --- Alpha --- %            
            if(~isempty(obj.R_Alpha))
                % Check if restricted parameters actually exist.
                if(r==0)
                    fprintf('\nWARNING: Cannot impose restrictions on alpha if r=0!\n')
                    fprintf('Either remove the restriction (set R_Alpha = [])  or increase rank (set r>0).\n'); 
                    error('Imposing restrictions on empty parameters');
                end
                % Check if column length of R_Alpha matches the number of
                % parameters.
                if(size(obj.R_Alpha,2) ~= p*r)
                    fprintf('\nWARNING: The length of R_Alpha does not match the number of parameters!\n');
                    fprintf('Please respecify R_Alpha so that the number of columns is p*r.\n');
                    error('Restriction misspecification');                  
                end
                % Check for redundancies.
                if(rank(obj.R_Alpha) < size(obj.R_Alpha,1))
                    fprintf('\nWARNING: R_Alpha has reduced rank!\n');
					fprintf('\nRedefine R_Alpha with linearly independent restrictions.\n');
                     error('Redundant restrictions in R_Alpha matrix');                  
                end
                % Check if user has imposed non-homogeneous alpha restrictions.
                if(~isempty(obj.r_Alpha) && (obj.r_Alpha ~= 0))
                    fprintf('\nWARNING: r_Alpha contains non-homogeneous restrictions (r_alpha non-zero).\n');
                    fprintf('All alpha restrictions have been made homogeneous.\n');
                end
                obj.r_Alpha = zeros(size(obj.R_Alpha,1),1);
            end
            
            % --- Beta --- %
            if(~isempty(obj.R_Beta))
                % Check if restricted parameters actually exist.
                if(r==0)
                    fprintf('\nWARNING: Cannot impose restrictions on Beta if r=0!\n')
                    fprintf('Either remove the restriction (set R_Beta = [])  or increase rank (set r>0).\n'); 
                    error('Imposing restrictions on empty parameters');
                end
                % Check if column length of R_Beta matches the number of
                % parameters (note p1 not p).
                if(size(obj.R_Beta,2) ~= p1*r)
                    fprintf('\nWARNING: The length of R_Beta does not match the number of parameters!\n');
                    fprintf('Please respecify R_Beta so that the number of columns is p1*r.\n');
                    error('Restriction misspecification');                  
                end
                % Check for redundancies.
                if(rank(obj.R_Beta) < size(obj.R_Beta,1))
                    fprintf('\nWARNING: R_Beta has reduced rank!\n');
                    fprintf('\nRedefine R_Beta with linearly independent restrictions.\n');
                    error('Redundant restrictions in R_Beta matrix');                  
                end

                % If non-zero restrictions haven't been specified by
                % the user, then set the r_Beta vector to zero's.
                if(isempty(obj.r_Beta))
                    obj.r_Beta = zeros(size(obj.R_Beta,1),1);
                else
                    % If user has specified restrictions, then check if
                    % the dimensions of LHS and RHS match. Note that if
                    % a restricted constant is being estimated, then it
                    % needs to be accounted for in the restrictions.
                    if(size(obj.r_Beta,1) ~= size(obj.R_Beta,1))
                        fprintf('\nWARNING: Row dimensions of R_Beta and r_Beta do not match!\n');
                        fprintf('Please redefine these matrices so that dimensions match.\n');
                        fprintf('Note: if a restricted constant has been included, \n');
                        fprintf('the Beta matrix has an additional row.\n');
                        error('All restrictions must be specified in the r_Beta variable');
                    end
                end
            end
            
            % Return the updated object of estimation options.
            newObj = obj; 
        end
        

    end
    
end

