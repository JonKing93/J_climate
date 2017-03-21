function[] = plotRuleN(s)
%% Plots the convergence (or lack thereof) of the Rule N significance tests
%
% -- In --
% s: A structure from an EOF_Analysis

% Plot the true confidence threshold
figure();
plot(s.iterTrueConf)
xlabel('Number of Monte Carlo Iterations');
ylabel('True Confidence Level of Eigenvalues');
title('True confidence level of iteratively selected significant eigenvalues.')

% Plot the significant eigenvalues for each iteration
figure()
plot(zscore(s.iterConfEigs));
title('Standardized Significance Threshold for Random Eigenvalues');
xlabel('Number of Monte Carlo Iterations');
ylabel('Standardized Significance');