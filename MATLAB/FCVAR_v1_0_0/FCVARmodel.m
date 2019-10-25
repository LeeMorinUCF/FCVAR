% Demonstrating the estimation of the fractionally cointegrated VAR model.
% Lee Morin
% August 21, 2011
%
% FCVARmodel demonstrates the estimation of the fractional cointegration model.
% 
% Calls the following functions:
% DefaultEstnOptions(db0, k, r, p) to set estimation options.
% RankTests(x, k, db0, FCVARoptions) to perform tests to determine
%       cointegrating rank.
% FCVARestn(x, k, r, FCVARoptions) to estimate the model without restrictions.
% GetResiduals(x, k, r, dbHat, alphaHat, betaHat, rhoHat, GammaHat, EstnOptions)
%       to compute and plot the residuals for diagnostics.
% RestrictEstnOptions(db0, k, r, p, defaultFCVARoptions) to reset estimation options
%       to impose various restrictions on the model.
% FCVARestn(x, k, r, db0, FCVARoptions) to repeat estimation with restrictions.
% GetResiduals(x, k, r, dbHat, alphaHat, betaHat, rhoHat, GammaHat, EstnOptions)
%       to compute and plot the residuals for further diagnostics.
% 
%______________________________________________________


clear

% Load data.
disp(' ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%  Loading Data  %%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
x = load('TheDanishData.csv');
p = size(x, 2); % Number of variables.

plot(x)
title('Data Used for Estimation')
pause




% Perform model selection beginning with lag length.
disp(' ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%  Model Selection: Size of db0  %%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

% Initial values.
% db0 = [ 1 ];        % One parameter implies restriction d = b.
db0 = [ 1 1 ];      % Two parameters define d and b separately.
if length(db0) == 1
    disp( 'Restriction d = b is imposed.');
elseif length(db0) == 2
    disp( 'Parameters d and b estimated separately.');
else
    disp( 'Warning: Starting values db0 improperly specified.');
end


% Perform model selection beginning with lag length.
disp(' ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%  Model Selection: Lag Length  %%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

% Select lag length (after performing a series of tests omitted here).
k = 2



% Now compute rank test stats.
disp(' ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%  Cointegrating Rank Tests  %%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')


% Obtain default options settings for unrestricted estimation.
rankTestOptions = DefaultEstnOptions(db0, k, p, p);

[ rankTestStats ] = rankTests(x, k, db0, rankTestOptions);

% Select cointegrating rank.
r = 2
pause


% Perform Unrestricted Estimation.
disp(' ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%  Unrestricted Estimation  %%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')


% Obtain default options settings for unrestricted estimation.
FCVARoptions = DefaultEstnOptions(db0, k, r, p);

% Perform estimation.
[ estimates ] = FCVARestn(x, k, r, db0, FCVARoptions);

% Obtain values from estimated parameters.
dbHat = cell2mat(estimates(1));
alphaHat = cell2mat(estimates(2));
betaHat = cell2mat(estimates(3));
rhoHat = cell2mat(estimates(4));
PiHat = cell2mat(estimates(5));
OmegaHat = cell2mat(estimates(6));
GammaHat = cell2mat(estimates(7));

pause


% Compute residuals from unrestricted estimation.
disp(' ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%  Computing Residuals  %%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

% Obtain EstnOptions from FCVARoptions.
EstnOptions = FCVARoptions{2};

% Now compute and plot residuals.
[ epsilon ] = GetResiduals(x, k, r, dbHat, alphaHat, betaHat, rhoHat, GammaHat, EstnOptions);
plot(epsilon);
title('Residuals from Estimation of Unrestricted Estimation.')
pause



% Perform Restricted Estimation.
disp(' ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%  Restricted Estimation  %%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

% Set optimization settings.

% Initial values.
% db0 = [ 1 ];        % One parameter implies restriction d = b.
db0 = [ 1 1 ];      % Two parameters define d and b separately.

% Obtain default options settings.
restictedFCVARoptions = RestrictEstnOptions(db0, k, r, p, FCVARoptions);

% Perform estimation.
[ restrictedEstimates ] = FCVARestn(x, k, r, restrictedFCVARoptions);

% Obtain values from estimated parameters.
dbHatR = cell2mat(restrictedEstimates(1));
alphaHatR = cell2mat(restrictedEstimates(2));
betaHatR = cell2mat(restrictedEstimates(3));
rhoHatR = cell2mat(restrictedEstimates(4));
PiHatR = cell2mat(restrictedEstimates(5));
OmegaHatR = cell2mat(restrictedEstimates(6));
GammaHatR = cell2mat(restrictedEstimates(7));

pause


% Compute residuals from Restricted estimation.
disp(' ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%  Computing Residuals  %%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

% Obtain EstnOptions from FCVARoptions.
restrictedEstnOptions = FCVARoptions{2};

% Now compute and plot residuals.
[ epsilonR ] = GetResiduals(x, k, r, dbHat, alphaHatR, betaHatR, rhoHatR, GammaHatR, restrictedEstnOptions);
plot(epsilon);
title('Residuals from Estimation of Restricted Estimation.')
pause



% The End.
disp(' ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%  The End  %%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')


% end
