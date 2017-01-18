function[SSA_Data] = SSA_Analysis(Data, M, algorithm, MC, noise, pval)
%% Performs singular spectra analyses for a collection of time series.
% Computes RCs (reconstructed components), and performs a Monte Carlo
% significance test
%
% -------------------------------------------------------------------------
% INPUTS:
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
% -------------------------------------------------------------------------
% OUTPUTS:
%
% [SSA_Data]: a structure with the following fields...THIS IS NOT
% FINISHED!!!
%
% eigvals: The eigenvalues of the time series. Each column corresponds to
%   a particular time series.
%
% eigvecs: The eigenvectors of each time series. Each dim1 x dim2 matrix
%   represents the vectors for a particular time series. Each column of
%   each matrix contains 1 eigenvector.
%
% RC: The reconstructed components of each time series. Each dim1 x dim2
%   matrix corresponds to a particular time series. Each column of each
%   matrix contains 1 reconstructed component
%
% sigEigDex: A boolean vector containing true/false indices for the
%   eigenvalues that passed the MC_SSA significance test.
%
% highSurrTail: A vector containing the surrogate eigenvalues on the upper
%   tail of the MC_SSA
%
% lowSurrTail: A vector containing the surrogate eigenvalues on the low
%   tail of the MC_SSA


%% Error check, get some startup values, set the vector of window sizes
setup(Data, M);


%% Remove mean from data
Data_m0 = detrend(Data, 'constant');


%% Get trajectories and covariance
[traj, C] = buildTrajandCov(Data_m0, M, algorithm);


%% Run an svd of each series
[eigvals, eigvecs] = quickSVD(C);


%% Get the RCs
[RC, ~] = getRCs(eigvecs, traj, algorithm);


%% Run an MC_EOF 
surrEig = MC_SSA(Data_m0, eigvecs, MC, noise, M, algorithm);


%% Do a significance test using the surrogate eigenvectors
[sigEigDex, highSurrTail, lowSurrTail] = sigTestMCSSA(pval, eigvals, surrEig);


%% Get the periods and frequencies of the data and surrogate eigenvalue tails


%% Assign data to a structure
SSA_Data = struct('Eigenvalues',eigvals,'Eigenvectors',eigvecs,'RC',RC,...
    'sigEigIndex',sigEigDex,'highSurrTail',highSurrTail,'lowSurrTail',lowSurrTail);
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