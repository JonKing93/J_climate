function[s] = SSA_Analysis(ts, M, algorithm, MC, noise, pval)
%% Performs singular spectra analyses for a collection of time series.
% Computes RCs (reconstructed components), and performs a Monte Carlo
% significance test
%
% [s] = SSA_Analysis(ts, M, algorithm, MC, noise, pval)
% Conducts a singular spectrum analysis of a time series and performs an
% MC-SSA significance test on the results.
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
%   surrEig: The surrogate eigenvalues from the Monte Carlo test. dim1 is
%       the Monte Carlo number, dim2 is the eigenvalue, dim3 is the time
%       series.
%
%   sigEigDex: A boolean vector containing true/false indices for the
%       eigenvalues that passed the MC_SSA significance test.
%
%   upperTail: A vector containing the surrogate eigenvalues from the
%       MC_SSA on the high end of the significance tail.
%
%   lowerTail: A vector containing the surrogate eigenvalues from the
%       MC_SSA on the low end of the significance tail
%       tail of the MC_SSA
%
%   maxFreq: The frequency with maximum power (using a periodogram) for
%       each Data eigenvector.
%
%   maxPeriod: The period with maximum power (using a periodogram) for each
%       Data eigenvector.
%
%   metadata: A cell with the trajectory algorithm, Monte Carlo number,
%       noise type, and p value for the analysis.

% Declare the initial structure
s = struct();

% Run an SSA on the time series
[s.singVals, s.singVecs, s.ts_m0, s.T, s.C] = simpleSSA(ts, M, algorithm);

% Get the RCs
[s.RCs, ~] = getRCs( s.singVecs, s.T, algorithm);

%% Run an MC_SSA
s.surrEig = MC_SSA(s.Data_m0, s.eigvecs, MC, noise, M, algorithm);

%% Do a significance test using the surrogate eigenvectors
[s.sigEigDex, s.upperTail, s.lowerTail] = sigTestMCSSA(pval, s.eigvals, s.surrEig);

%% Get the periods/frequencies with maximum power from the data eigenvectors
s.maxFreq = NaN(size(s.eigvals));
s.maxPeriod = s.maxFreq;
for k = 1:size(s.eigvals,2)
    [s.maxFreq, s.maxPeriod] = maxFreqPeriod( s.eigvecs(:,:,k));
end

%% Assign the metadata
s.metadata = {algorithm, MC, noise, pval};

end
