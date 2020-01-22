% This script demonstrates some of the additional features of the FCVAR
%   software package:
%
%   - Forecasting
%   - Bootstrap test of hypothesis on model coefficients
%   - Bootstrap rank test
%   - Simulation of fractionally cointegrated process
%
%  This script calls the estimates, option settings and data from
%   replication_JNP2014.m, so that file should be run before the examples
%   presented here.


%% --------- FORECAST ---------- %

% Forecast from the final restricted model.
NumPeriods = 12; % forecast horizon set to 12 months ahead.

% Assign the model whose coefficients will be used for forecasting.
modelF = m1r4;

xf = FCVARforecast(x1, modelF, NumPeriods);

model = m1r4;

% Series including forecast.
seriesF = [x1; xf];

% Equilibrium relation including forecasts.
equilF = seriesF*modelF.coeffs.betaHat;

% Determine the size of the vertical line to delimit data and forecast
%   values.
T = size(x1,1);
yMaxS  = max(max(seriesF));
yMinS  = min(min(seriesF));
yMaxEq = max(max(equilF));
yMinEq = min(min(equilF));

% Plot the results.
figure
subplot(2,1,1);
plot(seriesF),
title('Series including forecast'), xlabel('t');
line([T T], [yMinS yMaxS], 'Color','k');
subplot(2,1,2);
plot(equilF),
title('Equilibrium relation including forecasts'), xlabel('t');
line([T T], [yMinEq yMaxEq], 'Color','k');


%% --------- BOOTSTRAP HYPOTHESIS TEST ---------- %

% Test restriction that political variables do not enter the
%   cointegrating relation(s).

% Turn off plots for bootstrapping.
DefaultOpt.plotRoots = 0;

% Define estimation options for unrestricted model (alternative)
optUNR = DefaultOpt;

% Define estimation options for restricted model (null)
optRES = DefaultOpt;
optRES.R_Beta = [1 0 0];

% Number of bootstrap samples to generate
B = 999;

% Call to open the distributed processing (comment out if unavailable)
% matlabpool ('open',4); % for versions 2013a and earlier.
% parpool; % for versions 2013b and later.

[LRbs, H, mBS, mUNR] = FCVARboot(x1, k, r, optRES, optUNR, B);


%% Compare the bootstrap distribution to chi-squared distribution

% Estimate kernel density
[F,XI]=ksdensity(LRbs);

% Plot bootstrap density with chi-squared density
figure; plot(XI,F, XI, chi2pdf(XI,H.df))

legend(['Bootstrap PDF with ', num2str(B), ' BS samples'],...
    ['Chi Squared with ', num2str(H.df),' df'])


%% --------- BOOTSTRAP RANK TEST ---------- %

% Test rank 0 against rank 1
r1 = 0;
r2 = 1;

[LR_Rnk, H_Rnk, mBSr1, mBSr2] = FCVARbootRank(x1, k, DefaultOpt, r1, r2, B);

% Compare to P-value based on asymptotic distribution
fprintf('P-value: \t %1.3f\n', rankTestStats.pv(1));

% Close distributed processing (comment out if unavailable)
% matlabpool close;

%% --------- SIMULATION ---------- %

% Simulate the final restricted model, the same one used for forecasting
%   above.

% Number of periods to simulate
T_sim = 100;

% Simulate data
xSim = FCVARsim(x1, modelF, T_sim);

% Plot the results
figure;
plot(xSim)
legend('Support', 'Unemployment', 'Interest rate')
