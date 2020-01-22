function xf = FCVARforecast(x, model, NumPeriods)
% function xf = FCVARforecast(data, model, NumPeriods)
% Written by Michal Popiel and Morten Nielsen (This version 11.17.2014)
%
% DESCRIPTION: This function calculates recursive forecasts. It uses
%   FracDiff() and Lbk(), which are nested below.
%
% Input = data (T x p matrix of data)
%         model (a Matlab structure containing estimation results)
%         NumPeriods (number of steps ahead for forecast)
% Output = xf (NumPeriods x p matrix of forecasted values)
%_________________________________________________________________________

% --- Preliminary steps --- %
% x = data;
% p = size(data,2);
p = size(x,2);
opt = model.options;
cf  = model.coeffs;
d = cf.db(1);
b = cf.db(2);

% --- Recursively generate forecasts --- %
for i = 1:NumPeriods
    % append x with zeros to simplify calculations.
    x = [x; zeros(1,p)];
    T = size(x,1);

    % Adjust by level parameter if present.
    if(opt.levelParam)
        y = x - ones(T,1)*cf.muHat;
    else
        y = x;
    end

    % Main term, take fractional lag.
    z = Lbk(y,d,1);

    % Error correction term.
    if(~isempty(cf.alphaHat))
        z = z + FracDiff( Lbk(y, b, 1), d-b ) * cf.PiHat';
        if(opt.rConstant)
           z = z + FracDiff( Lbk(ones(T,1), b, 1), d-b ) * cf.rhoHat*cf.alphaHat';
        end
    end

    % Add unrestricted constant if present.
    if(opt.unrConstant)
        z = z + ones(T,1)*cf.xiHat';
    end
    
    % Add lags if present.
    if(~isempty(cf.GammaHat))
        k = size(cf.GammaHat,2)/p;
        z = z +  FracDiff(  Lbk( y , b, k)  , d) * cf.GammaHat';
    end

    % Adjust by level parameter if present.
    if(opt.levelParam)
        z = z + ones(T,1)*cf.muHat;
    end

    % Append forecast to x matrix.
    x = [x(1:T-1,:); z(T,:)];
end
    % Return forecasts.
    xf = x(size(data,1)+1:end,:);
end
