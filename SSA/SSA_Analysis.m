function[s] = SSA_Analysis(Data, M, algorithm, MC, noise, pval)
%% Performs singular spectra analyses for a collection of time series.
% Computes RCs (reconstructed components), and performs a Monte Carlo
% significance test
%
% [SSA] = SSA_Analysis(Data, M, algorithm, MC, noise, pval)
%
% ----- Inputs -----
%
% Data: a collection of time series stored as column vectors
%
% M: The size of the sampling window. If M is singular, all time series
% will be analyzed with a value of M. If a vector with length equal to the
% number of time series, the sampling window M_j will be applied to time
% series j.
%
% algorithm:
%   'BK': Broomhead-King
%   'VG': Vautard-Ghill
%
% MC: The monte carlo number for the MC_SSA significance test
%
% noise: The noise type to use in the MC_SSA significance test
%   'white': white Gaussian noise
%   'red': red, lag-1 autocorrelated noise
%
% pval: The desired confidence interval eigenvalues should test in the
%   MC_SSA significance test
%
% ----- Outputs -----
%
% s: a structure with the following fields
%
%   Data_m0: The time series normalized to a mean of 0
%
%   eigvals: The singular values of the time series. Each column 
%       corresponds to a particular time series.
%
%   eigvecs: The eigenvectors of each time series. Each dim1 x dim2 matrix
%       represents the vectors for a particular time series. Each column of
%       each matrix contains 1 eigenvector.
%
%   traj: The trajectory matrices constructed for the analyses. Each dim1 x
%       dim2 matrix is the trajectory matrix for one time series.
%
%   RCs: The reconstructed components of each time series. Each dim1 x dim2
%       matrix corresponds to a particular time series. Each column of each
%       matrix contains 1 reconstructed component.
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
%   metadata: A cell with (in this order) the trajectory algorithm, Monte
%       Carlo number, noise type, and p value for the analysis.


%% Error check, get some startup values, set the vector of window sizes
setup(Data, M);
s = struct();

%% Run a simple SSA
[s.eigvals, s.eigvecs, s.Data_m0, s.traj] = simpleSSA(Data, M, algorithm);

%% Get the RCs
[s.RCs, ~] = getRCs( s.eigvecs, s.traj, algorithm);

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

% ----- Helper functions -----
function[] = setup(Data, M)

% Ensure that Data is not 3D or greater
if ~ismatrix(Data)
    error('SSA_Analysis is for 2D collections of time series only.');
end

% Ensure M is singular
if length(M) ~= 1
    error('M must be singular');
end

end