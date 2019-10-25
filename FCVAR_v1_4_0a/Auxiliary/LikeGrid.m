function [ params ] = LikeGrid(x,k,r,opt)
% function [ params ] = LikeGrid(x,k,r,opt)
% Written by Michal Popiel and Morten Nielsen (This version 04.12.2016)
% 
% DESCRIPTION: This function evaluates the likelihood over a grid of values
% 	for (d,b) (or phi). It can be used when parameter estimates are sensitive to 
% 	starting values to give an approximation of the global max which can 
% 	then be used as the starting value in the numerical optimization in 
% 	FCVARestn().
%
% Input = x   (matrix of variables to be included in the system)
%         k   (number of lags)
%         r   (number of cointegrating vectors)
%         opt (object containing the estimation options)
% Output = params (row vector of d,b, and mu (if level parameter is selected)
%					corresponding to a maximum over the grid of (d,b), or phi)
%
% Note:	If opt.LocalMax == 0, LikeGrid returns the parameter values
%       corresponding to the global maximum of the likelihood on the grid.
%       If opt.LocalMax == 1, LikeGrid returns the parameter values for the
%       local maximum corresponding to the highest value of b. This
%       alleviates the identification problem mentioned in Johansen and
%       Nielsen (2010, section 2.3).
%_________________________________________________________________________

    p = size(x,2);
    
    if(opt.progress ~= 0)
        if(opt.progress == 1)
            M_status_bar = waitbar(0,['Model: k= ',num2str(k), ', r= ', num2str(r)]);
        end
        lastTic = tic();
    end
        
    
	
    % --- INITIALIZATION --- %
    % Check if performing a 2-dimensional or 1-dimensional grid search. The
    %   step size will be smaller if searching in only one dimension.
    if(isempty(opt.R_psi)) 
        Grid2d = 1;
        dbStep = 0.02;
    else
        Grid2d = 0;
        dbStep = 0.01;
    end
    
    % If equality restrictions are imposed, need to construct a grid over
    %   phi and obtain db = H*phi + h. 
    if(~isempty(opt.R_psi))
        % This set of variables is defined for easy translation between
        %   phi (unrestricted parameter) and d,b (restricted parameters).
        R = opt.R_psi;
        s = opt.r_psi;
        H = null(R);
        h = R'*inv(R*R')*s; 
    end      

    % GetBounds returns upper and lower bounds in terms of phi when
    % restrictions are imposed and in terms of d and b when they are
    % unrestricted. It also adjusts for limits on d and b imposed by user.
    [dbMax, dbMin] = GetBounds(opt);
    
    % Set up the grid.
    if(Grid2d)
        % Search over d as well as b.
        % Set up grid for d parameter
        dGrid =  dbMin(1):dbStep:dbMax(1);
        % Number of grid points
        nD    = length(dGrid);
        % Set up grid for b parameter
        bGrid =  dbMin(2):dbStep:dbMax(2);
        % Number of grid points
        nB    = length(bGrid);
        if(opt.constrained)
            % Since it is possible to have different grid lengths for d and
            % b when the parameters have different max/min values, the
            % following function calculates the total number of iterations
            % by counting the number of instances in which the grid values
            % in b are less than or equal to those in d. This should be
            % moved to the beginning of the loop for faster computation
            % later.
            totIters = sum(sum(bsxfun(@le,repmat(bGrid',1, nD),dGrid)));
        else
            % Unconstrained so search through entire grid for d.
            dStart   = 1;
            totIters = nB*nD;
        end
    else
        % Only searching over one parameter.
        bGrid =  dbMin(1):dbStep:dbMax(1);
        nB    = length(bGrid); 
        nD    = 1;
        dStart   = 1;
        totIters = nB;
    end


    % Create a matrix of NaN's to store likelihoods, we use NaN's here
    %   because NaN entries are not plotted and do not affect the finding
    %   of the maximum.
    like  = ones(nB,nD)*NaN;

    % Initialize storage bin and starting values for optimization involving 
    %   level parameter.
    if(opt.levelParam)
        mu  = zeros(p, nB, nD);
        StartVal = x(1,:);
    end

    iterCount = 0;
	
    % --- CALCULATE LIKELIHOOD OVER THE GRID --- %
    for iB=1:nB
        
        b = bGrid(iB);
        
        % If d>=b (constrained) then search only over values of d that are
        %   greater than or equal to b. Also, perform this operation only if
        %   searching over both d and b.
        if(opt.constrained && Grid2d)
            % Start at the index that corresponds to the first value in the
            % grid for d that is >= b
            dStart =  find(dGrid>=b, 1);
        end
        
        for iD=dStart:nD
           
            iterCount = iterCount + 1;
            

            % db is definied in terms of d,b for use by FCVARlikeMU
            % if level parameters are present and for displaying in
            % the output. phi is used by FCVARlike, which can
            % handle both phi or d,b and makes appropriate
            % adjustments inside the function.
            if(isempty(opt.R_psi))
                d = dGrid(iD);
                db = [ d b ];
                phi = db;
            else
                phi = bGrid(iB);
                db = H*phi + h;
            end
        
            if(opt.levelParam)
                % Optimize over level parameter for given (d,b).
                [ muHat, maxLike, ~ ] ...
                    = fminunc(@( params ) -FCVARlikeMu(x, db, params, k, r, opt), ...
                    StartVal, opt.UncFminOptions );
                % Store the results.
                like(iB,iD) = -maxLike;
                mu(:, iB,iD) = muHat;
            else
                % Only called if no level parameters. phi contains
                % either one or two parameters, depending on
                % whether or not restrictions have been imposed.
                like(iB,iD)   = FCVARlike(x, phi, k, r, opt);
            end
            
            if(opt.progress ~= 0)
                if(toc(lastTic) > opt.updateTime || iterCount == totIters) 
                    if(opt.progress == 1)
                        SimNotes = sprintf('Model: k=%g, r=%g\nb=%4.2f, d=%4.2f, like=%8.2f',...
                                    k, r, db(2),db(1), like(iB,iD) );                
                        waitbar(iterCount/totIters,M_status_bar, SimNotes);
                    else
                        fprintf('Progress = %5.1f%%, b=%4.2f, d=%4.2f, like=%g\n',...
                                    (iterCount/totIters)*100, db(2), db(1), like(iB,iD) );                
                    end
                    lastTic = tic();
                end          
            end
          
        end
   
    end


    % --- FIND THE MAX OVER THE GRID --- %

    if(opt.LocalMax)
    % Local max.
        if(Grid2d)
            [~,smax,~,~] = extrema2(like,1);
            [indexB,indexD] = ind2sub(size(like),smax);
        else
            [~,indexB,~,~] = extrema(like);
            indexD = 1;
        end
        % If there is no local max, return global max.
        if(isempty(indexB) || isempty(indexD))
            [ indexB, indexD ] = find(like == max(max(like)));
        end
    else
    % Global max.
        [ indexB, indexD ] = find(like == max(max(like)));
    end
    
    if(length(indexD)>1 || length(indexB)>1)
        % If maximum is not unique, take the index pair corresponding to 
        % the highest value of b.
        if(Grid2d)
            % Sort in ascending order according to indexB.
            [indexB,indBindx] = sort(indexB);
            indexD = indexD(indBindx);
        else
            % Sort in ascending order.
            indexB = sort(indexB);
        end
        indexB = indexB(end);
        indexD = indexD(end);
        fprintf('\nWarning, grid search did not find a unique maximum of the log-likelihood function.\n')
    end
          
    % Translate to d,b.
    if(~isempty(opt.R_psi))
        dbHatStar = (H*bGrid(indexB) + h)';
    else
        dbHatStar = [dGrid(indexD) bGrid(indexB)];
    end

    % --- STORE THE PARAMETER VALUES -- %
    params = dbHatStar;

    % Add level parameter corresponding to max likelihood.
    if(opt.levelParam)
        muHatStar  = mu(:, indexB, indexD)';
        params = [params muHatStar];
    end

    
    % --- PLOT THE LIKELIHOOD --- %
    if opt.plotLike
        figure;
        if(Grid2d)
            % 2-dimensional plot.
            mesh(dGrid, bGrid,like), xlabel('d'), ...
                zlabel('log-likelihood'), ylabel('b'), ...
                title(['Rank: ',num2str(r),' Lag: ',num2str(k)]);
        else
            % 1-dimensional plot.
            plot(bGrid, like) , 
            if(isempty(opt.R_psi))
                xlabel('d=b'),
            else
                xlabel('phi'), 
            end
            ...
                ylabel('log-likelihood'), ...
                title(['Rank: ',num2str(r),' Lag: ',num2str(k)]);
        end
    end   
end
