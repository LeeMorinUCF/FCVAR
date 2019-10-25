function [ UB, LB ] = GetBounds(opt)
% function [ UB, LB ] = GetBounds(opt)
% Written by Michal Popiel and Morten Nielsen (This version 04.11.2016)
% 
% DESCRIPTION: This function obtains upper and lower bounds on d,b or on 
%   phi, given by db = H*phi + h. 
%
% Input  = opt   (object containing estimation options)
% Output = UB (a 2x1 or 1x1 upper bound for db or phi)
%          LB (a 2x1 or 1x1 lower bound for db or phi)
%_________________________________________________________________________

if(isempty(opt.R_psi))
    % R_psi empty. Upper and lower bounds are the max and min values input
    % by the user. Both different and same bounds are permitted for d, b.
    % Maximum
    UB = opt.dbMax;
    % Minimum
    LB = opt.dbMin;
else
    % This set of variables is defined for easy translation between
    % phi (unrestricted parameter) and d,b (restricted parameters).
    R = opt.R_psi;
    s = opt.r_psi;
    H = null(R);
    h = R'*inv(R*R')*s; 

    UB = [];
    LB = [];

    % If there are two restrictions imposed, both d and b are restricted
    % and the upper and lower bounds are set to their restricted values.
    if(size(R,1) == 2)
        UB = h';
        LB = h';
    end

    % If there is 1 restriction, then there is 1 free parameter.
    if(size(R,1) == 1)        
        % Calculate endpoints of the grid from dbMin and dbMax to free
        % parameter phi. Note that since the null space can be either
        % positive or negative, the following conditional statements are
        % needed.
        if(H(1)>0)
            phiMin1 = (opt.dbMin(1) - h(1)) / H(1);
            phiMax1 = (opt.dbMax(1) - h(1)) / H(1);
        elseif(H(1)<0)
            phiMin1 = (opt.dbMax(1) - h(1)) / H(1);
            phiMax1 = (opt.dbMin(1) - h(1)) / H(1);
        else
            phiMin1 = NaN;
            phiMax1 = NaN;
        end
        if(H(2)>0)
            phiMin2 = (opt.dbMin(2) - h(2)) / H(2);
            phiMax2 = (opt.dbMax(2) - h(2)) / H(2);
        elseif(H(2)<0)
            phiMin2 = (opt.dbMax(2) - h(2)) / H(2);
            phiMax2 = (opt.dbMin(2) - h(2)) / H(2);
        else
            phiMin2 = NaN;
            phiMax2 = NaN;
        end

        % Take into account the condition d>=b, if required.
        if(opt.constrained)
            if (H(1)>H(2))
                phiMin3 = (h(2)-h(1))/(H(1) - H(2));
                phiMax3 = NaN;
            else
                phiMax3 = (h(2)-h(1))/(H(1) - H(2));
                phiMin3 = NaN;
            end
        else
            phiMin3 = NaN;
            phiMax3 = NaN;
        end

        LB = max(max(phiMin1, phiMin2), phiMin3);
        UB = min(min(phiMax1, phiMax2), phiMax3);    
    end
    
end   


end
