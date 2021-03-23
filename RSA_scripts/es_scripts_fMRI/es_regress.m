function [betas,R2] = es_regress(Y,X)
% Regression for multiple dependent variables
% Quicker than using Matlab's regress() in a for-loop
% Returns betas and R2 statistic
% Ed Sohoglu 2020

betas = pinv(X)*Y;
prediction = X*betas;
SS_effect = sum((prediction-mean(prediction)).^2); % effect sum of squares
SS_total = sum((Y-mean(Y)).^2); % total sum of squares
R2 = SS_effect./SS_total; % R2 statistic i.e. proportion of variance explained