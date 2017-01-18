function[eigvals, eigvecs, Data_m0, traj, C] = simpleSSA(Data, M, algorithm, nEigs)
%% Runs a single-spectrum analysis.
%
%--------------------------------------------------------------------------
% INPUTS:
%
% Data: A matrix of time series.
%   Each column of Data will be treated as an individual time series for
%   SSA
%
% M: The window size used to construct the trajectory matrix
%
% algorithm: the method used to construct the trajectory matrix
%     'BK':     (default) Broomhead-King
%     'VG':     Vautard-Ghil
%
% nEigs: (Optional Parameter) The number of eigenvalues and eigenvectors to
%   obtain. Will use a svds rather than svd analysis. (GENERALLY NOT
%   RECOMMENDED).
%
%--------------------------------------------------------------------------
% OUTPUTS:
%
% eigvals: The eigenvalues of each SSA of each time series
%
% eigvecs: The eigenvectors of each SSA of each time series 
%
% Data_m0: a set of time series with means removed
%
% traj: A set of trajectory matrices for the time series. Each dim1 x dim2
% matrix corresponds to a time series
%
% C: The covariance matrix for each trajectory matrix
%

%% Error check
setup(Data, M);


%% Remove mean from data
Data_m0 = detrend(Data, 'constant');


%% Get trajectories and covariance
[traj, C] = buildTrajandCov(Data_m0, M, algorithm);


%% Run an svd of each series
[eigvals, eigvecs] = quickSVD(C);


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