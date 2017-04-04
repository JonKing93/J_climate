function[s] = SSA_Analysis(ts, M, algorithm, MC, noise, pval, varargin)
%% Performs singular spectra analyses for a collection of time series.
% Computes RCs (reconstructed components), and performs a Monte Carlo
% significance test
%
% [s] = SSA_Analysis(ts, M, algorithm, MC, noise, pval)
% Conducts a singular spectrum analysis of a time series and performs an
% MC-SSA significance test on the results.
%
% [s] = SSA_Analysis(..., 'showProgress')
% A flag to display the current Monte Carlo iteration number in the
% significance test.
%
% [s] = SSA_Analysis(..., 'noConvergeTest')
% A flag to block convergence testing during the significance test. This
% may speed large analyses, but will cause a loss of information.
%
%
% ----- Inputs -----
%
% ts: A time series vector with equally spaced observations.
%
% M: The embedding dimension. This is the size of the sampling window used 
%       to construct trajectories in the data series. M corresponds to the 
%       largest wavelength that can be extracted from the data. However,
%       smaller values of M yield higher confidence in extracted information.
%
% algorithm:
%       'BK': Broomhead-King -- Slightly less bias for nonstationary time series
%       'VG': Vautard-Ghil -- Enhanced noise reduction for short time series
%
% MC: The number of iterations for the monte carlo significance test.
%
% noise: The noise type to use in the MC_SSA significance test
%       'white': white Gaussian noise
%       'red': red, lag-1 autocorrelated noise with added white noise
%
% pval: The desired significance level that eigenvalues should pass in the
%       MC_SSA significance test
%
%
% ----- Outputs -----
%
% s: a structure with the following fields
%
%   ts_m0: The time series normalized to a mean of 0
%
%   singVals: The singular values of the time series.
%
%   singVecs: The singular vectors of each time series. Each column is a
%       singular vector.
%
%   T: The trajectory matrix constructed for the analyses. Each dim1 x
%       dim2 matrix is the trajectory matrix for one time series.
%
%   C: The covariance matrix for the trajectory matrix.
%
%   expvar: The explained variance of each eigenvalue / eigenmode
%
%   RCs: The reconstructed components of each time series. Each dim1 x dim2
%       matrix corresponds to a particular time series. Each column of each
%       matrix contains 1 reconstructed component.
%
%   surrVals: The sorted surrogate eigenvalues from the Monte Carlo
%       significance test.
%
%   isSigVal: A boolean vector containing true/false indices for the
%       singular values that passed the MC_SSA significance test.
%
%   upSigVals: A vector containing the surrogate eigenvalues from the
%       MC_SSA on the upper tail of the confidence interval.
%
%   lowSigVals: A vector containing the surrogate eigenvalues from the
%       MC_SSA on the lower tail of the confidence interval.
%
%   maxFreq: The frequency with maximum power (using a raw periodogram) for
%       each singular vector.
%
%   maxPeriod: The period with maximum power (using a raw periodogram) for each
%       singular vector.
%
%   metadata: A cell with the trajectory algorithm, Monte Carlo number,
%       noise type, and p value for the analysis.

% Parse the inputs
[showProgress, convergeTest] = parseInputs(varargin{:});

% Declare the initial structure
s = struct();

% Run an SSA on the time series
[s.singVals, s.singVecs, s.ts_m0, s.T, s.C] = simpleSSA(ts, M, algorithm);

% Get the RCs
[s.RCs, ~] = getRCs( s.singVecs, s.T, algorithm);

% Run an MC_SSA
[s.surrVals, s.iterSigVals, s.iterTrueConf] = ...
    MC_SSA(s.ts_m0, s.singVecs, MC, noise, M, algorithm, pval, showProgress, convergeTest);

% Do a significance test using the surrogate eigenvectors
[s.iSigVal, s.upSigVals, s.lowSigVals] = sigTestMCSSA(pval, s.singVals, s.surrVals);

% Get the periods/frequencies with maximum power from the data eigenvectors
s.maxFreq = NaN(size(s.eigvals));
s.maxPeriod = s.maxFreq;
for k = 1:size(s.eigvals,2)
    [s.maxFreq, s.maxPeriod] = maxFreqPeriod( s.eigvecs(:,:,k));
end

%% Assign the metadata
s.metadata = {algorithm, MC, noise, pval};

end

function[showProgress, convergeTest] = parseInputs(varargin)
inArgs = varargin;
showProgress = 'noProgress';
convergeTest = 'convergeTest';

if ~isempty(inArgs)
    for k = 1:length(inArgs)
        arg = inArgs{k};
        if strcmpi(arg, 'showProgress')
            showProgress = arg;
        elseif strcmpi(arg, 'noConvergeTst')
            convergeTest = arg;
        else
            error('Unrecognized input');
        end
    end
end
end
