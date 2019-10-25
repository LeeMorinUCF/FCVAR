function [ FCVARoptions ] = AdjustEstnOptions(db0, k, r, p, initialFCVARoptions)
% 
% function [ FCVARoptions ] = AdjustEstnOptions(db0, k, r, p, initialFCVARoptions)
% Lee Morin & Morten Nielsen
% August 21, 2011
%
% 
% AdjustEstnOptions(db0, k, r, p, initialFCVARoptions) completes the specification 
%   of the options array, initialized by DefaultEstnOptions(db0, k, r, p)
%   or RestrictEstnOptions(db0, k, r, p, defaultFCVARoptions).
% 
% 
% input = vector or scalar db0 of fractional differencing parameters.
%       The size of db0 indicates whether or not restriction d = b is imposed
%           based on the number of free parameters.
%       scalars k, r and p denote the number of lags, cointegrating rank
%           and the number of variables, respectively.
%       array of arrays initialFCVARoptions contains the preliminary
%           specification of FCVARoptions initialized by DefaultEstnOptions(db0, k, r, p)
%           or RestrictEstnOptions(db0, k, r, p, defaultFCVARoptions).
% 
% 
% output = cell array FCVARoptions of options for use in estimation routines.
% 
% FCVARoptions = { dbEstnOptions, EstnOptions }, an array of 2 arrays of options
%       with the same form as initialFCVARoptions.
% 
% dbEstnOptions is an array of options for the estimation of the fractional differencing parameters d and b.
% dbEstnOptions = { constrained, C_db, c_db, R_db, r_db, dbFminOptions };
%   constrained is an indicator for the restriction that dbMax >= d >= b > dbMin.
%   C_db and c_db are the matrices that specify the above restriction.
%   dbMax is an upper bound for the estimation of d and b.
%   dbMin is a lower bound on d and b to convert a restriction of the form 
%       d >= b >= dbMin to represent the restriction d >= b > 0 with strict
%       inequality.
%   R_db and r_db are optional matrices for defining restrictions of the form R_db*db = r_db.
%   dbFminOptions is an array of options for built-in matlab functions.
% 
%   EstnOptions is an array of options for estimating the remaining parameters.
%   EstnOptions = { N, M, deterministics, R_Gamma, r_Gamma, [], [], [],
%                   nCols, A, H, a, b, h, [], [], [], [], 
%                   [], [], [], print2screen, printGammas,
%                   restrictFminOptions }
%       which is changed to 
%   EstnOptions = { N, M, deterministics, R_Gamma, r_Gamma, [], [], [],
%                   nCols, A, H, a, b, h, Abar, Aperp, replaceCols, keepCols,
%                   [], [], [], print2screen, printGammas,
%                   restrictFminOptions }
%
%   A and H are replaced with identities where no restrictions are specified.
%   Abar is an array of matrices such that Abar_i'*Abar = I, except in some special cases.
%   Aperp is an array of matrices such that Aperp_i*A_i = 0, except in some special cases.
%   replaceCols and keepCols are arrays of vectors specifying the column
%       numbers that are replaced or assumed known in each stage of the switching algorithm.
% 
%   h extends the restriction on beta to a similar restriction on
%       betaStar = (beta', rho)'.
% 
%   print2screen is a scalar indicator to print results (from RankTests and FCVARestn 
%      to screen and to generate graphs of characteristic roots.
%   printGammas is a scalar indicator to print coefficients of lagged differences as well.
% 
%   restrictFminOptions is an array of options for optimization over alpha and beta in the 
%       switching algorithm implemented in GetParamsSwitching.
%______________________________________________________



% Set optimization settings for the first stage optimization (d and b).
dbEstnOptions = initialFCVARoptions{1};

% Constrain optimization so that dbMax > d >= b > 0.
constrained = dbEstnOptions{1};
dbMax = dbEstnOptions{2};
dbMin = dbEstnOptions{3};

% For default case, dbMax should be the scalar dbMax, 
%   then the matrices C_db and c_db will be generated here
%   and will replace dbMax and dbMin,
%   else for restricted case, C_db and c_db are already adjusted.
if (numel(dbMax) == 1) && (numel(dbMin) == 1)
    
    if ( length(db0) == 1 ) && constrained
        C_db = [ -1; 1 ];
        c_db = [ -dbMin; dbMax ];
    elseif ( length(db0) == 2 ) && constrained
        C_db = [ -1 1; 0 -1; 1 0 ];
        c_db = [ 0; -dbMin; dbMax ];
    else
        C_db = [];
        c_db = [];
    end
    % Constrains d and b so that C_db*dbHat <= c_db.

    dbEstnOptions{2} = C_db;
    dbEstnOptions{3} = c_db;
end






% Now generate restriction matrices corresponding to primary specifications.
EstnOptions = initialFCVARoptions{2};
% EstnOptions = { N, M, deterministics, R_Gamma, r_Gamma, [], [], [],...
%                   nCols, A, H, a, b, h, Abar, Aperp, replaceCols, keepCols, 
%                   [], [], [], print2screen, printGammas, RestFminOptions };



% Set restrictions on Gamma (pk x p) matrix if required.

deterministics = EstnOptions{3};
R_Gamma = EstnOptions{4};
if ~isempty(R_Gamma)
    r_Gamma = zeros(size(R_Gamma, 1), 1);
    EstnOptions{5} = r_Gamma;
end
% Note: Does not allow user to estimate nonzero restrictions.




% Now adjust the array of arrays, { nCols, A, H, a, b, h }, each element of which has blanks
% wherever not specified.

% Recall:
% EstnOptions = { N, M, deterministics, R_Gamma, r_Gamma, [], [], [],...
%                   nCols, A, H, a, b, h, Abar, Aperp, replaceCols, keepCols, 
%                   [], [], [], print2screen, printGammas, RestFminOptions };
nCols = EstnOptions{9};
nRestns = length(nCols);    % Note: nRestns = 0 when nCols = [].
A = EstnOptions{10};
H = EstnOptions{11};
a = EstnOptions{12};
b = EstnOptions{13};
h = EstnOptions{14};
% To be filled in from blanks.
% Note: All are arrays same lengths as H above.
Abar = EstnOptions{15};
Aperp = EstnOptions{16};
replaceCols = EstnOptions{17};
keepCols = EstnOptions{18};


% Modify the specification of restrictions.
% Loop on j for usual case (special cases later).
% Note: nRestns = 0 when nCols = [] so this step is skipped if no restrictions are specified.
for j = 1:nRestns 
    
    % Hypotheses of the form alpha = (A_1*psi_1,...A_rA*psi_rA)
    %   where A_i is p x m_i, psi_i is m_i x r_i and sum(r_i) = r.
    if isempty(A(j))
        A{j} = eye(p);
    end
    
    % Hypotheses of the form beta = (H_1*phi_1,...H_rH*phi_rH)
    %   where H_i is p x s_i, phi_i is s_i x r_i and sum(r_i) = r.
    if isempty(H(j))
        H{j} = eye(p);
    end
    
    % Choose restrictions for known parameters.
    % a and b allow for a row of alpha or a column of beta to be fixed.
    % Note that known restrictions override linear restrictions in case the
    % user specifies both.
    
    
    % Calculate M-P inverse of Aj (unless different).
    Abar{j} = A{j}*inv(A{j}'*A{j});

    % Calculate matrix spanning orthogonal subspace of Aj.
    if A{j} == eye(p)
        % Blank if Aj = eye(p) (ie: no restrictions).
        Aperp{j} = [];
    else
        % AjPerp  = M_{Aj} = eye(p) - AjBar*Aj'.
        Aperp{j} = eye(p) - A{j}*inv(A{j}'*A{j})*A{j}'; 
    end
    
    % Determine relevant columns of beta0 to keep or drop on each iteration.
    replaceCols{j} = sum(nCols( (1:r) < j )) + 1 : sum(nCols( (1:r) < j + 1 ));
    keepCols{j} = [ ( 1 : sum(nCols( (1:r) < j )) )     ( sum(nCols( (1:r) < j + 1 )) + 1 : r ) ];
    
end






% Incorporate restrictions on rho into restrictions on betaStar.
% Two cases depend on whether or not each column of beta is known or has linear
% restrictions.
% In the case of linear restrictions, a vector of ones in h{j} indicates
% that rho{j} is to be estimated freely.
% All other cases impose the linear restriction using the same linear combination
% as for the corresponding columns of beta (ie: append to the columns of H{j}).
if strcmp(deterministics, 'restricted constant') 
    for j = 1:nRestns
        
        if isempty(b{j})
            % Append h's onto H's to create Hstar.
            % This works the same as H for beta in the switching algorithm.
            if h{j} == ones(size(h{j})) % rho estimated freely.
                H{j} = blkdiag( H{j}, diag(h{j}) );
            elseif h{j} == zeros(size(h{j})) % rho restricted to zero.
                H{j} = [ H{j}; h{j} ];
            else % rho restricted to same linear comination as beta.
                H{j} = [ H{j}; h{j} ];
            end
        else
            % Append h's onto b's to restrict entire column of betaStar.
            b{j} = [ b{j}; h{j} ]; % rho specified as known along with beta.
        end
        
    end
end






% Consider special cases in Johansen's book: no iterations required.
% For all other cases, switching algorithm is primed as usual.
% Here it is understood that the single known set of columns is represented
%       in the first set, otherwise, the switching algorithm is used as is.

% Special cases in Johansen's book.
if     length(nCols) == 2 &&  (  ~isempty(a(1)) &&  isempty(a(2))  )...  
         && (   A(1) == eye(p) && A(2) == eye(p) ...
            &&  B(2) == eye(p) && B(2) == eye(p) ...
            &&  isempty(b(1))  &&  isempty(b(2))    )
    
    % H(2) already set to identities, while Aperp(2) already blank.
    % Algorithm looks for alpha in space orthogonal to a(1).
    % A(2) set to a(1)perp.
    A(2) = a{1}*inv(a{1}'*a{1});
    % Abar(2) set to (a(1)perp)Bar.
    Abar(2) = eye(p) - a{1}*inv(a{1}'*a{1})*a{1}';
	
    
elseif length(nCols) == 2 &&  (  ~isempty(b(1)) &&  isempty(b(2))  )...
         && (   A(1) == eye(p) && A(2) == eye(p) ...
            &&  B(2) == eye(p) && B(2) == eye(p) ...
            &&  isempty(a(1))  &&  isempty(a(2))    )
    
    % H(2), A(2) and Abar(2) already set to identities.
    % Algorithm looks for beta in space orthogonal to b(1).
    Aperp(2) = b(1);
    
end    % Default :Keep settings as above.






% Now replace and create the restriction arrays as outputs.
EstnOptions{10} = A;
EstnOptions{11} = H;
% Unchanged:
% a = EstnOptions{12};
% b = EstnOptions{13};
% h = EstnOptions{14};
% Newly created:
EstnOptions{15} = Abar;
EstnOptions{16} = Aperp;
EstnOptions{17} = replaceCols;
EstnOptions{18} = keepCols;




% Now put together the adjusted version of the options array.
% dbEstnOptions carries the options for the first stage estimation of fractional parameters.
% EstnOptions carries the options for all other functions and restrictions.
FCVARoptions = { dbEstnOptions, EstnOptions};



% end