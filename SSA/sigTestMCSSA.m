function[isSigVal, upSigVals, lowSigVals] = sigTestMCSSA( pval, singVals, surrVals)
%% Performs a significance test of data SSA eigenvalues against MC_SSA surrogates eigenvalues
%
% [sigEigDex, upperTail, lowerTail] = sigTestMCSSA( pval, dataEigvals, surrEigvals)
% Determines which singular values pass the significance threshold.
%
% ----- Inputs -----
%
% pval: The significance level that the test should pass. This should be a
%   number between 0 and 1.
%
% dsingVals: A set of singular values.
%
% surrVals: A set of surrogate eigenvalues.
%
%
% ----- Outputs -----
%
% isSigVal: A boolean (true/false) vector with the indices of the data 
%   eigenvalues that pass the significance test. Each column corresponds to a
%   particular time series.
%
% upSigVals: A matrix of the surrogate eigenvalues on the upper tail of
%   the significance test. Each column corresponds to a particular time
%   series.
%
% lowSigVals: A matrix of the surrogate eigenvalues on the lower tail of
%   the significance test. Each column corresponds to a particular time
%   series.
%

% Run an error check and get some initial sizes.
[MC] = setup(pval, singVals, surrVals);
conf = 1-pval;

% Get the indices of the surrogate band.
maxRowDex = ceil( MC * (conf + (pval/2)  ) );
minRowDex = floor( MC * (pval/2) );

% Get the vectors of the two tails of the confidence interval;
upSigVals = surrVals( maxRowDex, :);
lowSigVals = surrVals( minRowDex, :);

% Get the indices of the eigenvalues that passed the significance test
isSigVal = singVals > upSigVals;

end


% ----- Helper functions -----
function[MC] = setup(pval, singVals, surrVals)

% Ensure pval is between 0 and 1
if pval >= 1 || pval <= 0
    error('p Value must be between 0 and 1');
end

% Ensure the inputs have the correct dimensionality
if ~isvector(singVals)
    error('singVals must be a vector');
end

if ~ismatrix(surrVals)
    error('surrVals must be a matrix');
end

% Ensure the dimensions of the inputs align
[dwindow] = length(singVals);
[MC, swindow] = size(surrVals);

if dwindow ~= swindow
    error('Data and Surrogate eigenvalue matrices are not correctly aligned');
end

% Ensure the desired confidence interval is possible for the MC number
if floor( MC * pval/2 ) == 0
    minMC = 1/(pval/2);
    error('The Monte Carlo number is too small for this confidence interval. A value of %i or larger is required.',minMC);
end
end